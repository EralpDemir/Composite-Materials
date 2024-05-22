%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is written by Eralp Demir                                   %
% to design the optimized stiffness of the fiber-steered                  %
% variable-stiffness composite materials                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%close all
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
[Asen, Bsen, Dsen] = OPT_tangent(t, U);

%__________________________ OPTIMIZATION ________________________________
% Only symmetric laminate has been considered

% Initial guess
All_xi_A0 = zeros(4,numel);
All_xi_D0 = zeros(4,numel);

% % Use these initial conditions for in-plane cases!
% All_xi_A0(1,:) = 1;
% All_xi_A0(3,:) = 1;

% % Use these initial conditions for out-of-plane cases!
% All_xi_D0(2,:) = 0.25;
% All_xi_D0(4,:) = 0.25;

% Symmetric type of laminates are considered only!
All_xi_B0 = zeros(4,numel);
All_xi_B = All_xi_B0;
All_xi_A=zeros(4,numel);
All_xi_D=zeros(4,numel);


% Run FEM analysis for the initial values of Lamination Parameters
[d, e, r] = FEM_lam(numel, nnpe, dof, numind, nqpts, dofst, ndir, ...
    conn, crds, All_xi_A0, All_xi_B0, All_xi_D0, t, U, dofBCs, ...
    valBCs, numBCs, f, K, qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, qL_u, ...
    Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp);

disp(['initial strain energy:  '  num2str(e)])




% Constant stiffness optimization
if strcmp(typSTIFF, 'CON')
    
    % Set up the optimization variable
    switch typOPT
       

        case 1 % in-plane problem
            
            
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
                    
                    


                      
                              
                    x0 = zeros(4,1);

                    A = [ 0 0  1 0;
                          0 0 -1 0 ];

                    b = [ 1;
                          1;];
                      
                      

            end

        case 2 % out-of-plane problem
            
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
                

                  
                  
                x0 = zeros(4,1);

                
                
                A = [ 0 0  1 0;
                      0 0 -1 0 ];

                b = [ 1;
                      1;];

                

                  
        end

    end

    % equality and lower-upeer bounds of the optimization problem
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];

    % Non-linear constraint function
    nonlcon = @(y)OPT_constraint_con(y,typOPT,typLAM);
    
    objective = @(x)OPT_objective_con( x, numel, nnpe, dof, numind, nqpts, ...
        dofst, ndir, conn, crds, All_xi_A0, All_xi_B0, All_xi_D0, t, U, ...
        dofBCs, valBCs, numBCs, f, K, qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, ...
        qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp, Asen, Bsen, ...
        Dsen, typOPT, typLAM);
    
    
    
    [x, funcval, ef, output, alpha] = ...
        fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    
    % Assign the optimized results
    if typOPT == 1
        
        switch typLAM
            
            case 'ORT'
                                
                for i=1:1:numel
                    All_xi_A(1,i) = x(1);
                    All_xi_A(3,i) = x(2);
                end
                All_xi_D = All_xi_D0;
                
            case 'GEN'

                for i=1:1:numel
                    All_xi_A(1:4,i) = x;
                end
                All_xi_D = All_xi_D0;
                
        end

    elseif typOPT == 2
 
        switch typLAM
            
            case 'ORT'
                
                All_xi_A = All_xi_A0;
                for i=1:1:numel
                    All_xi_D(1,i) = x(1);
                    All_xi_D(3,i) = x(2);
                end
                
                
            case 'GEN'
                
                All_xi_A = All_xi_A0;
                for i=1:1:numel
                    All_xi_D(1:4,i) = x;
                end
                
        end

    end
    
    
    
% Varible stiffness optimization
elseif strcmp(typSTIFF, 'VAR')



    
    % Set up the optimization variable
    switch typOPT

        case 1 % in-plane problem
            
            
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
            
            
            

                    x0 = reshape(All_xi_A0, 4*numel,1);


                    A = zeros(2*numel,2*numel);
                    
                    b = zeros(2*numel,1);
                    
                    for iele = 1:1:numel
                        
                        A(2*iele-1:1:2*iele, 4*iele-3:1:4*iele) = [  0 0 1 0;
                            0 0 -1 0];

                        b(2*iele-1:2*iele, 1) = [ 1;
                            1;];

                    end
                    
                    
                    
            end

        case 2 % out-of-plane problem

            
            
            
            switch typLAM
                
                case 'ORT'
            
                    x0 = reshape(All_xi_D0(1:2:end,:), 2*numel,1);
                    
                    
                    
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
            
                    x0 = reshape(All_xi_D0, 4*numel,1);
                    
                    
                    
                    
                    
                    A = zeros(2*numel,4*numel);
                    
                    b = zeros(2*numel,1);
                    
                    
                    for iele = 1:1:numel


                                                
                                                
                        A(2*iele-1:1:2*iele, 4*iele-3:1:4*iele) = [  0 0 1 0;
                            0 0 -1 0];

                        b(2*iele-1:2*iele, 1) = [ 1;
                            1;];                        
                    end

                    
                
            end

    end

    % equality and lower-upeer bounds of the optimization problem
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];

    % Non-linear constraint function
    nonlcon = @(y)OPT_constraint_var(y,typOPT,numel,typLAM);

    objective = @(x)OPT_objective_var( x, numel, nnpe, dof, numind, nqpts, ...
        dofst, ndir, conn, crds, All_xi_A0, All_xi_B0, All_xi_D0, t, U, ...
        dofBCs, valBCs, numBCs, f, K, qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, ...
        qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp, Asen, Bsen, ...
        Dsen, typOPT, typLAM);

    [x, funcval, ef, output, alpha] = ...
        fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    % Assign the optimized results
    if typOPT == 1
        
        switch typLAM
            
            case 'ORT'
                                
                for i=1:1:numel
                    All_xi_A(1,i) = x(2*i-1);
                    All_xi_A(3,i) = x(2*i);
                end
                All_xi_D = All_xi_D0;
                
            case 'GEN'

                All_xi_A = reshape(x,4,numel);
                All_xi_D = All_xi_D0;
        

        end
                
    elseif typOPT == 2
        
         switch typLAM
            
            case 'ORT'
                
                All_xi_A = All_xi_A0;
                for i=1:1:numel
                    All_xi_D(1,i) = x(2*i-1);
                    All_xi_D(3,i) = x(2*i);
                end
                
                
            case 'GEN'

                All_xi_A = All_xi_A0;
                All_xi_D = reshape(x,4,numel);
                
         end
        
    end

end
    
    
All_angPly = zeros(numel,nPlys);
All_xi_A_lib = zeros(4,numel);
All_xi_D_lib = zeros(4,numel);

% Find the nearest lamination parameter in the space
for iele = 1:1:numel
    
    xi = [All_xi_A(1:4,iele); All_xi_D(1:4,iele);];
    
    [ang, xi_A, xi_D] = LP_find(xi_As, xi_Ds, angPlys, noAngPoss, xi, typOPT);
    
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

