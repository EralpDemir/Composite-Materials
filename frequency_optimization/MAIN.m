%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is written by Eralp Demir, and Ali Rashed to design         %
% the optimized natural frequency of the fiber-steered                    %
% variable-stiffness composite materials                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

warning('off') %#ok<WNOFF>
format long

%__________________________ INITIALIZATIONS ______________________________
% FEM Inputs
FEM_inputs;

% Inputs material properties
MAT_inputs;

% Inputs for laminate specifications
LP_inputs;

% Inputs for the optimization procedure
OPT_inputs;

% FEM computation for least square fit and continuity constraints
LSC_inputs;

% Initialize LSTCC constant
LSC_initialize;

% Construct the Lamination Parameter Space
[xi_As, xi_Bs, xi_Ds, noAngPoss, angPlys] = ...
    LP_space(ang_range, t, Ltype, nPlys, zPly);

% Calculate analytical 1st order sensitivities
[Asen, Dsen] = OPT_tangent(t, U);

%__________________________ OPTIMIZATION ________________________________
% Only symmetric laminate has been considered

% Initial guess
All_xi_A0 = zeros(4,numel);
All_xi_D0 = zeros(4,numel);

%[All_xi_A0, All_xi_D0] = random_inputs(numel);







% Initialize lamination parameters
All_xi_A = All_xi_A0;
All_xi_D = All_xi_D0;



% Run FEM analysis for the initial values of Lamination Parameters
[lambda0, a0] = FEM_lam(numel, nnpe, dof, conn, crds, rho, ...
    All_xi_A, All_xi_D, t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp);

% Natural frequency non-dimensional
disp(['initial frequency:  '  num2str(sqrt(lambda0) * nondim)])




% Constant stiffness optimization
if strcmp(typSTIFF, 'CON')     

    
    
    
     switch typLAM
                
                case 'ORT'
                    

                   x0 = zeros(2,1);

                    A = [  1  0;
                          -1  0;
                           0  1;
                           0 -1];

                    b = [ 1;
                          1;
                          1;
                          1];
            
                                  
                      
                      
                case 'GEN'
    
    
            
                    x0 = [  All_xi_D0(1,1);
                            All_xi_D0(2,1);
                            All_xi_D0(3,1);
                            All_xi_D0(4,1);];



                    A = [ 0 0  1 0;
                          0 0 -1 0 ];

                    b = [ 1;
                          1;];


     end




    % equality and lower-upeer bounds of the optimization problem
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];

    % Non-linear constraint function
    nonlcon = @(y)OPT_constraint_con(y,typLAM);
    
    objective = @(x)OPT_objective_con(x, numel, nnpe, dof, conn, crds, rho, ...
        t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp, Dsen, typLAM);
    
    
    
    
    
    
    [x, funcval, ef, output, alpha] = ...
        fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    
    % Assign the optimized results

        
    switch typLAM

        case 'ORT'

            for i=1:1:numel
                All_xi_D(1,i) = x(1);
                All_xi_D(3,i) = x(2);
            end
            All_xi_A = All_xi_A0;

        case 'GEN'

            for i=1:1:numel
                All_xi_D(1:4,i) = x;
            end
            All_xi_A = All_xi_A0;

    end



    


    
    
% Varible stiffness optimization
elseif strcmp(typSTIFF, 'VAR')

    switch typLAM
                
        case 'ORT'



            x0 = reshape(All_xi_A0(1:2:end,:), 2*numel,1);



            A = zeros(4*numel, 2*numel);

            b = zeros(4*numel, 1);

            for iele = 1:1:numel





                A(4*iele-3:1:4*iele, 2*iele-1:1:2*iele) = [  1  0;
                                                            -1  0;
                                                             0  1;
                                                             0 -1; ];

                b(4*iele-3:4*iele, 1) = [   1;
                                            1;
                                            1;
                                            1;];

            end





        case 'GEN'

    
            % Set up the optimization variable
            x0 = reshape(All_xi_D0, 4*numel,1);


            A = zeros(2*numel, 4*numel);

            b = zeros(2*numel, 1);

            for iele = 1:1:numel

                                A(2*iele-1:1:2*iele, 4*iele-3:1:4*iele) = [ 0 0 1 0;
                                                                            0 0 -1 0];

                                b(2*iele-1:2*iele, 1) = [   1;
                                                            1;];

            end

    end
     

    % equality and lower-upper bounds of the optimization problem
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];


    % Non-linear constraint function
    nonlcon = @(y)OPT_constraint_var(y,numel,typLAM);

    objective = @(x)OPT_objective_var(x, numel, nnpe, dof, conn, crds, rho, ...
        t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp, Dsen, typLAM);

    [x, funcval] = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    % Assign the optimized results
    switch typLAM
            
            case 'ORT'
                                
                for i=1:1:numel
                    All_xi_D(1,i) = x(2*i-1);
                    All_xi_D(3,i) = x(2*i);
                end
                All_xi_A = All_xi_A0;
                
            case 'GEN'

                All_xi_A = zeros(4,numel);
                All_xi_D = zeros(4,numel);

                for i=1:1:numel
                    All_xi_D(1,i) = x(4*i-3,1);
                    All_xi_D(2,i) = x(4*i-2,1);
                    All_xi_D(3,i) = x(4*i-1,1);
                    All_xi_D(4,i) = x(4*i,1);
                end
                
    end

 

end
   


% Run FEM analysis for the final/optimized values of Lamination Parameters
[lambda, a] = FEM_lam(numel, nnpe, dof, conn, crds, rho, ...
    All_xi_A, All_xi_D, t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp);



% Natural frequency non-dimensional
disp(['optimized frequency:  '  num2str(sqrt(lambda) * nondim)])


All_angPly = zeros(numel,nPlys);
All_xi_A_lib = zeros(4,numel);
All_xi_D_lib = zeros(4,numel);

% Find the nearest lamination parameter in the space
for iele = 1:1:numel
    
    xi = [All_xi_D(1,iele); All_xi_D(2,iele); All_xi_D(3,iele); All_xi_D(4,iele);];
    
    [ang, xi_A, xi_D] = LP_find(xi_As, xi_Ds, angPlys, noAngPoss, xi);
    
    All_angPly(iele,1:nPlys) = ang;
    
    All_xi_A_lib(1:4,iele) = xi_A;
    
    All_xi_D_lib(1:4,iele) = xi_D;
    
end

% Smooth the angles
[angPlyNodes] = LSC_angdist(BigNmat, Bigfvec, numel, LSC_nnpe, ...
    LSC_conn, LSC_totnod, nPlys, All_angPly);

% Calculate the angles at the element centers
[All_angPly_c] = LSC_angel0(numel, LSC_nnpe, LSC_conn, nPlys, ...
    angPlyNodes, N0);

% Compute overall curvature
kappa = zeros(1,nPlys);
for iPly = 1:1:nPlys
    angs(:,1) = angPlyNodes(:,iPly)*pi/180;
    kappa(iPly) = sqrt(angs'*BigCmat*angs/Area)/2;
end

clear angs

% Compute curvature elementally
kappa_el = zeros(numel,nPlys);
for iele = 1:1:numel
    for iPly = 1:1:nPlys
        nodes = LSC_conn(iele,:);
        angs = angPlyNodes(nodes,iPly)*pi/180;
        Cmat(:,:) = Cmat_els(iele,:,:);
        kappa_el(iele,iPly) = sqrt(angs'*Cmat*angs)/2;
    end
end

clear angs nodes iele iPly

