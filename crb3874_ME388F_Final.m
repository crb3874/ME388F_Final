% - 2D, 2G Eigenvalue Solver using Deterministic Transport (MOCs) - %
% - Reuse MOCs code from HW4 as a solver.
        
ps = Input(true);
tests05x = MoCEigenSolver( ps );

function ps = Input( flag )

    %%% - Eigensolver Inputs - %%%
    ps = struct();
    ps.tol = 10^(-6);
    ps.chi = [1,0];
    ps.maxit = 200;
    
    if flag
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Problem Specifications for MoC %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % - Group 1 (Fast) - %
        mocps1 = struct();

        mocps1.dxVector = 1*ones(25,1);
        mocps1.dyVector = 1*ones(25,1);

        % Reflector outer region
        mocps1.matArray = 2 * ones(25,25);

        % Fueled inner region
        mocps1.matArray(5:20,5:20) = 1;

        mocps1.nPolr = 16;
        mocps1.nAzim = 16;
        mocps1.raySpace = 0.1;

        mocps1.maxit              = 10000;
        mocps1.tol                = 10^(-6);

        mocps1.boundaryConditions = [ [ "Vacuum" "Vacuum" ];
                                  [ "Vacuum" "Vacuum" ] ];

        % - Group 2 (Thermal) - %                 
        mocps2 = struct();

        mocps2.dxVector = 1*ones(25,1);
        mocps2.dyVector = 1*ones(25,1);

        % Reflector outer region
        mocps2.matArray = 4 * ones(25,25);

        % Fueled inner region
        mocps2.matArray(5:20,5:20) = 3;

        mocps2.nPolr = 16;
        mocps2.nAzim = 16;
        mocps2.raySpace = 0.1;

        mocps2.maxit              = 10000;
        mocps2.tol                = 10^(-6);

        mocps2.boundaryConditions = [ [ "Vacuum" "Vacuum" ];
                                  [ "Vacuum" "Vacuum" ] ]; 

        %%% - END MOC INPUTS - %%%
        moc1 = MoC(mocps1);     % Generate mesh and ray data for
        moc2 = MoC(mocps2);     % each energy group.
        ps.MoCData = {moc1, moc2};
    end
    

end
