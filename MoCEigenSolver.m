classdef MoCEigenSolver
    
    properties
        
        MoCData         % - Cell array containing all MoC objects (each group).
        nx              % - Problem dimension (x)
        ny              % - Problem dimension (y)
        ng              % - Number of groups.
        k               % - Eigenvalue.
        chi             % - Fission spectrum.
        scatterKernel0  % - Arrays for group-to-group scattering (zeroth moment).
        scatterKernel1  % - Arrays for group-to-group scattering (first moment).
        scatterKernel2  % - Arrays for group-to-group scattering (second moment).
        fissSource      % - Fission source (same array size as scalar flux).
        outertol        % - Convergence criteria for eigensolve.
        innertol        % - Convergence criteria for inner iteration.
        maxouterit      % - Maximum inner (sweep) iterations.
        maxinnerit      % - Maximum power iterations.
        khist           % - Eigenvalue iteration history.
        
    end
    %--End Properties Block--%
    
    methods
        
        function ins = MoCEigenSolver( ps )
            
            % Instantiate - store properties from ps, create first source
            % update, run eigensolve routine.
            
            ins.MoCData = ps.MoCData;
            ins.outertol = ps.tol;
            ins.innertol = ps.MoCData{1}.ps.tol;
            ins.chi = ps.chi;
            ins.nx  = ins.MoCData{1}.mD.nx;
            ins.ny  = ins.MoCData{2}.mD.ny;
            ins.ng  = length(ins.MoCData);
            ins.maxinnerit = ins.MoCData{1}.ps.maxit;
            ins.maxouterit = ps.maxit;
            
            ins.scatterKernel0 = cell(ins.ng);
            ins.scatterKernel1 = cell(ins.ng);
            ins.scatterKernel2 = cell(ins.ng);
            
            % Generate scattering kernels.
            for i = 1:ins.ng
                
                ins.scatterKernel0{i} = struct();
                ins.scatterKernel1{i} = struct();
                
                mat = ins.MoCData{i}.ps.matArray;
                                         
                scatterData = arrayfun( @(m) ins.ScatterData(m) , mat );
                
                % Scattering - zeroth moment, within group
                s0in = { scatterData.s0in };
                ins.scatterKernel0{i}.s0in = cell2mat( reshape( s0in, ins.nx , ins.ny ) );
                
                % Scattering - zeroth moment, out of group.
                s0out = { scatterData.s0out };
                ins.scatterKernel0{i}.s0out = cell2mat( reshape( s0out, ins.nx , ins.ny ) );
                
                % Scattering - first moment, within group.
                s1in = { scatterData.s1in };
                ins.scatterKernel1{i}.s1in = cell2mat( reshape( s1in, ins.nx , ins.ny ) );
                
                % Scattering - first moment, out of group.
                s1out = { scatterData.s1out };
                ins.scatterKernel1{i}.s1out = cell2mat( reshape( s1out, ins.nx , ins.ny ) );
                
                % Scattering - second moment, within group.
                s2in = { scatterData.s2in };
                ins.scatterKernel2{i}.s2in = cell2mat( reshape( s2in, ins.nx , ins.ny ) );
                
                % Scattering - second moment, out of group.
                s2out = { scatterData.s2out };
                ins.scatterKernel2{i}.s2out = cell2mat( reshape( s2out, ins.nx , ins.ny ) );
                
            end
            
            % Initial guess.
            ins.k = 1;
            
            ins = ins.EigenSolve();
            
        end
        
        function ins = EigenSolve( ins )
            % Use power iteration (outer loop).
            
            % Initial fission source (initial scalar flux assumed to be
            % all ones).
            
            ins.fissSource = zeros(ins.nx, ins.ny);
            ins.fissSource = ins.fissSource + ins.MoCData{1}.mD.nusf;
            ins.fissSource = ins.fissSource + ins.MoCData{2}.mD.nusf;
            
            ins.khist = [];
            it = 0;
            
            while it < ins.maxouterit
                
                it = it+1;
                
                % Sweep fixed fission source until covergence.
                
                ins.MoCData = ins.MGMoCSolve(); % Iterates scatter source to converge this
                                                % fission source.
                
                % Compute new fission source.
                newfissSource = zeros(ins.nx, ins.ny);
                newfissSource = newfissSource + ins.MoCData{1}.mD.nusf .* ins.MoCData{1}.scalarFlux;
                newfissSource = newfissSource + ins.MoCData{2}.mD.nusf .* ins.MoCData{2}.scalarFlux;
                
                % Update k.
                newk = ins.k * norm(newfissSource(:),2)/norm(ins.fissSource(:),2);
                
                % Check for convergence.
                e = abs(newk-ins.k)/ins.k;
                
                ins.k = newk;
                ins.fissSource = newfissSource;
                ins.khist = [ins.khist, ins.k];
                
                if e < ins.outertol
                    break
                end
                
            end
            
            if e > ins.outertol
                fprintf("MoC Eigensolve failed to converge within " + int2str(ins.maxouterit) + " power iterations.")
                fprintf("Current eigenvalue: " + ins.k + ". Current residual: " + e);
            end
            
        end
        
        function MoCData = MGMoCSolve( ins )
            
            % Sweep a fixed fission source to convergence with multi-group
            % scattering - can largely reuse the MoC solver with some work
            % at the end of the loop for updating multigroup scatter
            % sources.
            
            % Initialize scalar flux -
            ins.MoCData{1}.scalarFlux = zeros(ins.nx , ins.ny);
            ins.MoCData{2}.scalarFlux = zeros(ins.nx , ins.ny);
            % moment11 = zeros(ins.nx , ins.ny);
            % moment12 = zeros(ins.nx , ins.ny);
            % moment21 = zeros(ins.nx , ins.ny);
            % moment22 = zeros(ins.nx , ins.ny);
            
            ins.MoCData{1}.source = ins.chi(1)*ins.fissSource/(4*pi*ins.k);
            ins.MoCData{2}.source = ins.chi(2)*ins.fissSource/(4*pi*ins.k);
            
            % Sweep rays over fixed fission source until convergence.
            it = 0;
            
            while it < ins.maxinnerit
                
                it = it+1;
                scalarFlux1 = ins.MoCData{1}.scalarFlux;
                scalarFlux2 = ins.MoCData{2}.scalarFlux;
                
                % No need for inner upscattering iterations here.
                
                % - FAST GROUP - %
                % Fast Group scatter source.
                ssourceFast = ins.scatterKernel0{1}.s0in .* ins.MoCData{1}.scalarFlux;
                ssourceFast = ssourceFast + ins.scatterKernel0{2}.s0out .* ins.MoCData{2}.scalarFlux;                
                
                % Update MoC object with new fission + scatter source.
                ins.MoCData{1}.source = (1/4/pi)*(ins.chi(1) * ins.fissSource/ins.k + ssourceFast);
                % ins.MoCData{1}.linsource =  (1/4/pi) * (ins.scatterKernel1{1}.s1in .* moment11 + ins.scatterKernel1{2}.s1out .* moment21);
                % ins.MoCData{1}.quadsource = (1/4/pi) * (ins.scatterKernel2{1}.s2in .* moment12 + ins.scatterKernel2{2}.s2out .* moment22);
                
                % Sweep rays for fast group.
                [ins.MoCData{1}.scalarFlux, ~, ~] = ins.MoCData{1}.SweepRays(); 
                
                % - THERMAL GROUP - %
                % Thermal Group scatter source.
                ssourceThermal = ins.scatterKernel0{2}.s0in .* ins.MoCData{2}.scalarFlux;
                ssourceThermal = ssourceThermal + ins.scatterKernel0{1}.s0out .* ins.MoCData{1}.scalarFlux;
                
                % Update MoC object new fission + scatter source.
                ins.MoCData{2}.source = (1/4/pi)*(ins.chi(2) * ins.fissSource/ins.k + ssourceThermal);
                % ins.MoCData{2}.linsource =  (1/4/pi) * (ins.scatterKernel1{2}.s1in .* moment21 + ins.scatterKernel1{1}.s1out .* moment11);
                % ins.MoCData{2}.quadsource = (1/4/pi) * (ins.scatterKernel2{2}.s2in .* moment22 + ins.scatterKernel2{1}.s2out .* moment12);
                
                % Sweep rays for thermal group.
                [ins.MoCData{2}.scalarFlux, ~, ~] = ins.MoCData{2}.SweepRays();
                
                % - EVALUATE INNER ITERATION CONVERGENCE - %
                dPhi            = scalarFlux1 - ins.MoCData{1}.scalarFlux;
                dPhi            = dPhi + (scalarFlux2 - ins.MoCData{2}.scalarFlux);
                
                e = norm(dPhi(:), 2);
                
                if e < ins.innertol
                    break
                end
                
            end
            
            if e >ins.innertol
                fprintf( "MoC failed to converge to the specified tolerance (" + ins.ps.tol + ") in " + int2str(ins.ps.maxit) + " iterations" );
                fprintf( "Current residual: " + e );
            end
            
            MoCData = ins.MoCData;

        end
        
        function ScatterData = ScatterData( ~ , matID )
            
            FuelData1 = struct( "s0in",  3.40438E-01 , "s0out" ,  1.63429E-02 , "s1in", 4.38057E-01 , "s1out" , 2.06865E-02 , "s2in" , 3.22608E-01 , "s2out" , -1.06319E-02);
            WaterData1 = struct("s0in",  3.29611E-01 , "s0out" ,  2.50280E-02 , "s1in", 5.96010E-01 , "s1out" , 3.23976E-02 , "s2in" , 4.24315E-01 , "s2out" , -1.65064E-02);
            FuelData2 = struct( "s0in",  8.37316E-01 , "s0out" ,  0.000000000 , "s1in", 7.72600E-01 , "s1out" , 0.000000000 , "s2in" , 3.87959E-01 , "s2out" , 0.0000000000);
            WaterData2 = struct("s0in",  1.05166E+00 , "s0out" ,  0.000000000 , "s1in", 1.10727E+00 , "s1out" , 0.000000000 , "s2in" , 5.58949E-01 , "s2out" , 0.0000000000);
            
            XS_Directory = { FuelData1 , WaterData1 , FuelData2 , WaterData2 };

            ScatterData = XS_Directory{ matID };
    
        end
  
        function FluxPlotter( ins , g )
            % Call FluxPlotter methods from group g.
            
            ins.MoCData{g}.FluxPlotter();
            
        end
        
    end
    %--End Methods Block--%
    
end