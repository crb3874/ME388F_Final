classdef MoC
    properties
        
        ps              
        mD              % meshData object
        rayBin          % Nested cell array contain rayData objects
        source          % Scalar flux source -> divied up to each ordinate
        linsource       % Linearly anisotropic source.
        quadsource      % Quadratically anisotropic source.
        scalarFlux      % Our solution
        moment1         % linear moment.
        moment2         % quadratic moment.
        avgOuter        % Avg scalar flux in outer region
        mocData         % General outputs of interest (iterations, etc.)
        
    end
    
    methods
        
        function ins        = MoC( ps )
            
            % Store inputs
            ins.ps = ps;
            
            % Generate meshData from ps.
            ins.mD = meshData( ps );
            
            % Seed rays
            ins.rayBin = ins.SeedRays(  );
            
        end
        
        function MoCSolve( ins )
            
            % Initialize fixed source -
            ins.source = ins.mD.is0;
            
            % Initialize scalar flux -
            ins.scalarFlux = zeros(ins.mD.nx , ins.mD.ny);
            
            % Sweep rays over fixed source until convergence.
            it = 0;
            
            while it < ins.ps.maxit
                
                it = it+1;
                scalarFlux0     = ins.scalarFlux;
                ins.scalarFlux  = ins.SweepRays(); %Sweep rays :)
                dPhi            = scalarFlux0 - ins.scalarFlux;
                
                e = norm(dPhi(:), 2);
                
                % Convergence check
                if e < ins.ps.tol
                    break
                end
                
                % Update source
                ins.source      = ins.SourceUpdate();
                
            end
            
            if e >ins.ps.tol
                fprintf("MoC failed to converge to the specified tolerance (" + ins.ps.tol + ") in " + int2str(ins.ps.maxit) + " iterations");
            end
            
            Outer = (ins.ps.matArray~=3);
            ins.avgOuter = mean(ins.scalarFlux(Outer));
            
        end
        
        function rayBin     = SeedRays( ins ) 
            
            % rayBin will store ray data for all azimuthal angles -->
            % rayBin{aziAng}{ray1,ray2,ray3...} [All rayData objects]
            
            rayBin = { {} }; % Initialize nested cell array
            
            % Grab stuff from mD;
            aziAng = ins.mD.dOrd.aziAng;
            x      = ins.mD.x;
            y      = ins.mD.y;
            ox     = ins.mD.dOrd.ox;
            oy     = ins.mD.dOrd.oy;
            oz     = ins.mD.dOrd.oz;
            del    = ins.mD.dOrd.raySpace;
            
            % Now - - - need to generate rayData objects.
            % For each azimuthal angle, generate a new column in
            % rayBin. Step along and seed rays according to the 
            % spacing specified in mD; Seed right then seed up.
   
            
            % Loop through azimuthal angles.
            for o = 1:length(aziAng)
                
                % We need to seed the rays with the following info:
                % xs , ys , ox , oy , ozV , trI0 , trJ0 , mD
                % i.e. -> starting point, direction, and INDEx
                % of initial material (as well as meshData). Make sure we 
                % keep track of these things in the following loops
                
                % Track rays... 
                r = 0;
                
                %Seed right going rays from the left, seed 
                %left going rays from the right.
                
                if ox(o) > 0
                    xs=x(1);
                    xbound = x(end);
                    idx = 1;
                elseif ox(o) < 0
                    xs=x(end);
                    xbound = 0;
                    idx = length(x)-1;
                end
                
                ys = y(1);
                jdx = 1;
            
                % Seed along x...
                while sign(ox(o)) * (xbound-xs) > 0
                    r = r+1;
                    [~,ix] = min( abs( xs-x ) ); % Find which cell we're in
                    if x(ix) > xs || ix == length(x) % Seeded from right boundary, a little awk
                        idx = ix-1;
                    else
                        idx = ix;
                    end

                    %Seed ray.
                    newRay = rayData( xs , ys , ox(o) , oy(o) , oz , idx , jdx , ins.mD );
                    rayBin{o}{r} = newRay;
                    
                    %Step in x
                    xs = xs + sign(ox(o))*abs(del/oy(o));
                    
                end
                
                % Reset x...
                if ox(o) > 0
                    xs=x(1);
                    idx = 1;
                elseif ox(o) < 0
                    xs=x(end);
                    idx = length(x)-1;
                end

                % Seed along y...
                while ys < y(end)
                    r = r+1;
                    
                   	[~,jx] = min( abs( ys-y ) ); % Find which cell we're in
                    if y(jx) > ys
                        jdx = jx-1;
                    else
                        jdx = jx;
                    end
                    
                    %Seed ray.
                    newRay = rayData( xs , ys , ox(o) , oy(o) , oz , idx , jdx , ins.mD );
                    rayBin{o}{r} = newRay;
                    
                    %Step in y
                    ys = ys + abs(del/ox(o));      
                    
                end
                
            end
            
        end
        
        function [scalarFlux , moment1 , moment2] = SweepRays( ins )
            
            %Initialize zero array for the scalar flux.
            scalarFlux = zeros(ins.mD.nx,ins.mD.ny);
            moment1 = zeros(ins.mD.nx,ins.mD.ny);
            moment2 = zeros(ins.mD.nx,ins.mD.ny);
            
            %Loop through angles --
            for o = 1:length(ins.rayBin)
                
                %Loop through rays --
                for r = 1:length(ins.rayBin{o})
                    
                    %Initialize track incoming flux from boundary.
                    psi0 = zeros(ins.ps.nPolr,1);
                              % For now, only handling vacuum boundaries.
                              % Very challenging to handle reflecting
                              % as rays that would reflect may not line up
                              % At least, we've done nothing to force that to
                              % be possible...
                    
                    ray = ins.rayBin{o}{r};
                    
                    %Loop through tracks --> Forward
                    for k = 1:+1:length(ray.v_trL)
                        
                        %Call in values from ray- cell index and track length
                        i = ray.v_trI(k);
                        j = ray.v_trJ(k);
                        d = ray.v_trL(k);
                        
                        %Grab source value and st for this cell.
                        q = ins.source(i,j);
                        % q = q + ins.linsource(i,j) * ins.mD.dOrd.aziAng(o);
                        % q = q + ins.quadsource(i,j) * (1/2)*(3*ins.mD.dOrd.aziAng(o)^2-1);
                        
                        st = ins.mD.st(i,j);
                        
                        %Calculate exiting flux
                        psi1 = q/st + (psi0 - q/st).*ray.v_trA(:,k);
                        
                        %Calculate cell-averaged flux
                        psiAvg = q/st + (psi0 - psi1)./ray.v_trTau(:,k);
                        
                        %psi1 and psiAvg now contain rows for each polar
                        %angle - we can now sum over each with the proper
                        %weight to update the scalar flux in cell (i,j);
                        
                        rayA     = ins.mD.dOrd.raySpace * d;
                        wei      = ins.mD.dOrd.aziWei(o)*ins.mD.dOrd.polWei;
                        
                        scalarFlux(i,j) = scalarFlux(i,j) + rayA*wei'*psiAvg;
                        moment1(i,j)    = moment1(i,j) + rayA*wei'*psiAvg*ins.mD.dOrd.aziAng(o);
                        moment2(i,j)    = moment2(i,j) + (0.5)*rayA*wei'*psiAvg*(3*ins.mD.dOrd.aziAng(o)^2-1);
                        
                        %Outgoing flux is incoming for next k.
                        psi0 = psi1;
                        
                    end
                    %-- End forward track loop

                    %Initialize track incoming flux from boundary.
                    psi0 = zeros(ins.ps.nPolr,1);
                    
                    %Loop through tracks <-- Backwards
                    for k = length(ray.v_trL):-1:1
                        
                        %Call in values from ray- cell index and track length
                        i = ray.v_trI(k);
                        j = ray.v_trJ(k);
                        d = ray.v_trL(k);
                        
                        %Grab source value and st for this cell.
                        q = ins.source(i,j);
                        % q = q + ins.linsource(i,j) * ins.mD.dOrd.aziAng(o);
                        % q = q + ins.quadsource(i,j) * (1/2)*(3*ins.mD.dOrd.aziAng(o)^2-1);
                        
                        st = ins.mD.st(i,j);
                        
                        %Calculate exiting flux
                        psi1 = q/st + (psi0 - q/st).*ray.v_trA(:,k);
                        
                        %Calculate cell-averaged flux
                        psiAvg = q/st + (psi0 - psi1)./ray.v_trTau(:,k);
                        
                        %psi1 and psiAvg now contain rows for each polar
                        %angle - we can now sum over each with the proper
                        %weight to update the scalar flux in cell (i,j);
                        
                        rayA     = ins.mD.dOrd.raySpace * d;
                        wei       = ins.mD.dOrd.aziWei(o)*ins.mD.dOrd.polWei;
                        
                        scalarFlux(i,j) = scalarFlux(i,j) + rayA*wei'*psiAvg;
                        moment1(i,j)    = moment1(i,j) + rayA*wei'*psiAvg*ins.mD.dOrd.aziAng(o);
                        moment2(i,j)    = moment2(i,j) + (0.5)*rayA*wei'*psiAvg*(3*ins.mD.dOrd.aziAng(o)^2-1);
                        
                        %Outgoing flux is incoming for next k.
                        psi0 = psi1;
                        
                    end
                    %-- End backward track loop
                    
                end
                %-- End ray loop
                
            end
            %-- End angle loop
            
        end
        
        function source     = SourceUpdate( ins )
            
            % Updates source using most recent approximation for the scalar
            % flux. 
            
            source = ins.mD.is0;
            source = source + (1/4/pi) * ins.mD.ss0 .* ins.scalarFlux;
            
            % Could build on this for fission/anisotropic scattering, etc.
            
        end
        
        function RayPlotter( ins )
            
            for o = 1:ins.ps.nAzim/2
                
                figure(o);
                % Plot grid --

                hold on
                
                [X,Y] = meshgrid(ins.mD.x, ins.mD.y);
                plot(X,Y,'k')
                plot(Y,X,'k')

                for i = 1:length(ins.rayBin{o})

                    h=plot(ins.rayBin{o}{i}.xint , ins.rayBin{o}{i}.yint , '-o');
                    title("Rays for azim " + int2str(o))
                    drawnow
                    pause(0.03)

                end
                
                xlabel("X (cm)")
                ylabel("Y (cm)")
                
                hold off;
                
            end
        end
        
        function FluxPlotter( ins )
            
            xaxis = (ins.mD.x(2:end) + ins.mD.x(1:end-1))/2;
            yaxis = (ins.mD.y(2:end) + ins.mD.y(1:end-1))/2;
            xyaxis = sqrt(xaxis.^2 + yaxis.^2);
            phidiag = diag(ins.scalarFlux);
            
            figure(97);
            %Full diagonal plot
            semilogy(xyaxis, phidiag);
            title("Flux along full diagonal")
            xlabel("XY [cm]")
            ylabel("\phi (n/cm^2s)")
            
            figure(98);
            %Half diagonal plot
            semilogy(xyaxis(1:end/2), phidiag(1:end/2));
            title("Flux along half diagonal")
            xlabel("XY [cm]")
            ylabel("\phi (n/cm^2s)")
            
            figure(99);
            %Surface plot
            
            surf(xaxis, yaxis, ins.scalarFlux);
            xlabel("X [cm]")
            ylabel("Y [cm]")
            zlabel("\phi (n/cm^2s)")
            
        end
        
        %--End Methods Block--%
    end

end