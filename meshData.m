classdef meshData
    properties
        
        x
        y
        dx
        dy
        nx
        ny
        st
        sf
        nusf
        ss0
        ss1
        is0
        dOrd
        
        %%%%%%%%% END PROPERTIES BLOCK %%%%%%%%%
        
    end
    
    methods
        
        function ins = meshData( ps ) %Constructor for instance 'ins'
            
            %%%%%%%%%%%%%%%% MESH/COORDINATE DATA %%%%%%%%%%%%%%%%%%%%%%%%%
            ins.dx = ps.dxVector;
            ins.dy = ps.dyVector;
            
            ins.nx = length(ins.dx);
            ins.ny = length(ins.dy);
            
            ins.x  = [0; cumsum(ins.dx)];
            ins.y  = [0; cumsum(ins.dy)];
            
            %%%%%%%%%%%%%%%% MATERIAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            crossSections = arrayfun( @(m) ins.CrossSections(m) , ps.matArray );

            % Total cross section
            st  = { crossSections.st };
            ins.st  = cell2mat( reshape( st , ins.nx , ins.ny ) ); 

            % Scattering cross section (zeroth moment, within group)
            ss0 = { crossSections.ss0 };
            ins.ss0 = cell2mat( reshape( ss0, ins.nx , ins.ny ) );
            
            % Scattering cross section (first moment)
            ss1 = { crossSections.ss1 };
            ins.ss1 = cell2mat( reshape( ss1, ins.nx , ins.ny ) );
            
            % Fission cross section
            sf = { crossSections.sf };
            ins.sf = cell2mat( reshape( sf  , ins.nx , ins.ny ) );
            
            % nu-fission data
            nusf = { crossSections.nusf };
            ins.nusf = cell2mat( reshape( nusf, ins.nx , ins.ny ) );
            
            %%%%%%%%%%%%%%% ANGULAR DISCRETIZATION %%%%%%%%%%%%%%%%%%%%%%%%
            ins.dOrd = ins.AngularDiscretize( ps.nPolr , ps.nAzim );
            ins.dOrd.raySpace = ps.raySpace;
            
        end
        
        function xsData = CrossSections( ~ , matID )
           
            %% Material IDs: [Fuel , Reflector , Scatterer , Pure Absorber , Absorber , Isotropic , Anisotropic , Air]
            %%               [   1 ,         2 ,         3 ,             4 ,        5 ,         6 ,           7 ,   8]
            
            FuelData1  = struct("st" , 3.60902E-01 , "sf" , 1.29920E-03 , "ss0" , 4.12091E-03 , "ss1" , 0.1  , "sa" , 4.16193E-03 , "nusf" , 3.38910E-03);
            WaterData1 = struct("st" , 3.54704E-01 , "sf" , 0.000000000 , "ss0" , 6.51894E-05 , "ss1" , 0.1  , "sa" , 6.51820E-05 , "nusf" , 0.000000000);
            FuelData2  = struct("st" , 8.77427E-01 , "sf" , 2.47341E-02 , "ss0" , 4.01106E-02 , "ss1" , 0.1  , "sa" , 4.00994E-02 , "nusf" , 6.02693E-02);
            WaterData2 = struct("st" , 1.05606E+00 , "sf" , 0.000000000 , "ss0" , 4.40145E-03 , "ss1" , 0.1  , "sa" , 4.40145E-03 , "nusf" , 0.000000000);

            XS_Directory = {FuelData1, WaterData1, FuelData2, WaterData2};

            xsData = XS_Directory{ matID };
    
        end
        
        function dOrd = AngularDiscretize( ~ , nPolr , nAzim )
    
            % Generate polar (legendre) discrete ordinates.

            if nPolr == 0

                % Discrete angles.

                polAng = [-1 ; 1];
                polWei = [1 ; 1];

            else

                % Generate nPol discrete ordinates using Legendre points/weights
                % for z-axis ordinates.

                syms o

                legendreN = legendreP( nPolr , o );

                % Quadrature points.
                lsgqPointsSym           = root( legendreN , o );
                polAng                  = double ( vpa( lsgqPointsSym ) );

                % Weights require the derivative.
                legendreNDeriv(o)       = diff( legendreN ) ; 
                legendreNDerivEval      = double( vpa( legendreNDeriv( polAng ) ) );

                % Compute weights.
                polWei                  = 2./(1-polAng.^2)./(legendreNDerivEval.^2);

            end

           % Generate azimuthal (equally spaced) discrete ordinates.

           aziAng = (2*(0:(nAzim-1)/2)+1) * pi/nAzim;
           aziWei = (2*pi/nAzim) * ones(nAzim,1);

           % Output data

           dOrd = struct;
           dOrd.polAng = polAng;
           dOrd.polWei = polWei;
           dOrd.aziAng = aziAng;
           dOrd.aziWei = aziWei;
           dOrd.ox = cos(aziAng);
           dOrd.oy = sin(aziAng);
           dOrd.oz = cos(polAng);

        end
        
        %%%%%%%%% END METHODS BLOCK %%%%%%%%%
        
    end
end