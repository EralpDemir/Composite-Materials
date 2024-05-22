%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is written by Eralp Demir to design                         %
% the optimized buckling of the fiber-steered                             %
% variable-stiffness composite materials                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all



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


% [All_xi_A0, All_xi_D0] = random_inputs(numel);
% 



% Initialize lamination parameters
All_xi_A = All_xi_A0;
All_xi_D = All_xi_D0;



% Run FEM analysis for the initial values of Lamination Parameters
[R_m_bottom0, lambda0, a0, d0] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
    conn, crds, All_xi_A0, All_xi_D0, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, ...
    dofBCs_m_bottom, dofBCs_m_top);




disp(['initial eigenfrequency:  '  num2str(lambda0)])




% Constant stiffness optimization
if strcmp(typSTIFF, 'CON')     

            
            

    x0=zeros(4,1);
    % Take the first elements values since it is constant
    x0(1) = All_xi_A0(1,1);
    x0(2) = All_xi_A0(3,1);
    x0(3) = All_xi_D0(1,1);
    x0(4) = All_xi_D0(3,1);

    

    % inequality constraints

    A = [ 1  0  0  0;
         -1  0  0  0;
         
          0  1  0  0;
          0 -1  0  0;

          0  0  1  0;
          0  0 -1  0;
          
          0  0  0  1;
          0  0  0 -1;];
      

    b = [ 1;
          1;
          1;
          1;
          1;
          1;
          1;
          1;];                






    % equality and lower-upeer bounds of the optimization problem
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];

    % Non-linear constraint function
    nonlcon = @(x)OPT_constraint_con(x);
    
    objective = @(x)OPT_objective_con(x, numel, nnpe, dof_b, dof_m, ...
    conn, crds, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, f_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, numind_m, ...
    numind_b, Asen, Dsen,dofBCs_m_bottom, dofBCs_m_top);



    
    [x, funcval, ef, output, alpha] = ...
        fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    
    % Assign the optimized results
    for i=1:1:numel
        All_xi_A(1,i) = x(1);
        All_xi_A(3,i) = x(2);
        All_xi_D(1,i) = x(3);
        All_xi_D(3,i) = x(4);
    end
            


    


    
    
% Varible stiffness optimization
elseif strcmp(typSTIFF, 'VAR')



    
    % Set up the optimization variable




    x0=zeros(4*numel,1);
    for i=1:1:numel
        x0(4*i-3,1) = All_xi_A0(1,i);
        x0(4*i-2,1) = All_xi_A0(3,i);
        x0(4*i-1,1) = All_xi_D0(1,i);
        x0(4*i,1)   = All_xi_D0(3,i);
    end




    

    A = zeros(8*numel, 4*numel);

    b = zeros(8*numel, 1);

    for iele = 1:1:numel

        A(8*iele-7:1:8*iele, 4*iele-3:1:4*iele) = [ 1 0 0 0;
                                                   -1 0 0 0;
                                                    0  1 0 0;
                                                    0 -1 0 0;
                                                    0 0  1 0;
                                                    0 0 -1 0;
                                                    0 0 0  1;
                                                    0 0 0 -1;];

        b(8*iele-7:8*iele, 1) = [ 1;
            1;
            1;
            1;
            1;
            1;
            1;
            1;];
    end





     

    % equality and lower-upper bounds of the optimization problem
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];

    % Non-linear constraint function
    nonlcon = @(x)OPT_constraint_var(x,numel);

    objective = @(x)OPT_objective_var(x, numel, nnpe, dof_b, dof_m, ...
    conn, crds, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, f_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, numind_m, ...
    numind_b, Asen, Dsen, BigNmat, Bigfvec, LSC_nnpe, ...
    LSC_conn, LSC_totnod, N0,dofBCs_m_bottom, dofBCs_m_top);



    [x, funcval] = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    % Assign the optimized results
    All_xi_A = zeros(4,numel);
    All_xi_D = zeros(4,numel);

    for i=1:1:numel
        All_xi_A(1,i) = x(4*i-3,1);
        All_xi_A(3,i) = x(4*i-2,1);
        All_xi_D(1,i) = x(4*i-1,1);
        All_xi_D(3,i) = x(4*i,1);
    end

 

end
   


% Run FEM analysis for the final/optimized values of Lamination Parameters
[R_m_bottom, lambda, a] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
    conn, crds, All_xi_A, All_xi_D, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, ...
    dofBCs_m_bottom, dofBCs_m_top);




    
All_angPly = zeros(numel,nPlys);
All_xi_A_lib = zeros(4,numel);
All_xi_D_lib = zeros(4,numel);

% Find the nearest lamination parameter in the space
for iele = 1:1:numel
    
    xi = [All_xi_A(1,iele); All_xi_A(3,iele); All_xi_D(1,iele); All_xi_D(3,iele);];

%    xi = zeros(4,1);
    
    [ang, xi_A, xi_D] = LP_find(xi_As, xi_Ds, angPlys, noAngPoss, xi);
    
    All_angPly(iele,1:nPlys) = ang;
    
%     [xi_A,xi_B,xi_D] = MAT_lp(8,[90, 0, 0, 90, 90, 0, 0, 90], z, t);
    
    All_xi_A_lib(1:4,iele) = xi_A;
    
    All_xi_D_lib(1:4,iele) = xi_D;
    

    
    
    
end




% Run FEM analysis for the final/optimized values of Lamination Parameters
[R_m_bottom_lib,lambda_lib, a_lib] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
    conn, crds, All_xi_A_lib, All_xi_D_lib, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, dofBCs_m_bottom, dofBCs_m_top);







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


