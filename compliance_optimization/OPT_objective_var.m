function [obj, grad_obj] = OPT_objective_var(x,numel, nnpe, dof, numind, ...
   nqpts, dofst, ndir, conn,crds, ...
    All_xi_A, All_xi_B, All_xi_D, t, U, dofBCs, valBCs, numBCs, f, K, ...
    qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, qL_u, Fipedlu_ElNod, qL_v, ...
    Fipedlv_ElNod, eltyp, Asen, Bsen, Dsen, typOPT, typLAM)

% Assign Stiffness parameters
B(1:6) = zeros(6,1);

% Set up the optimization variable
switch typLAM
    
    case 'ORT'
           
         if typOPT == 1

            All_xi_A(1,1:numel) = x(1:2:end);
            All_xi_A(3,1:numel) = x(2:2:end);

        elseif typOPT == 2

            All_xi_D(1,1:numel) = x(1:2:end);
            All_xi_D(3,1:numel) = x(2:2:end);

        end
        
            
    case 'GEN'
        
        if typOPT == 1

            All_xi_A = reshape(x,4,numel);

        elseif typOPT == 2

            All_xi_D = reshape(x,4,numel);

        end
        
        
end




% Run FEM analysis for the new values of Lamination Parameters
[d, obj] = FEM_lam(numel, nnpe, dof, numind, nqpts, dofst, ndir, conn,...
    crds, All_xi_A, All_xi_B, All_xi_D, t, U, dofBCs, valBCs, numBCs, ...
    f, K, qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, qL_u, Fipedlu_ElNod, qL_v, ...
    Fipedlv_ElNod, eltyp);

grad_obj = [];

for iele = 1:1:numel
    
    % Elemental displacements
    d_el = [];
    for j = 1:1:nnpe
        d_el = [d_el; d(dof*conn(iele,j)-dof+1:1:dof*conn(iele,j))];
    end
    
    % Compute the sensitivities:
    % For each lamination parameter (A, B, D, H)
    de_dxi = zeros(4,4);
    for k = 1:2:4
        
        % For each lamination parameter type (1,2,3,4)
        for j = 1:1:4
            
            % Assign Stiffness parameters
            A = zeros(6,1);
            % B = zeros(6,1);
            D = zeros(6,1);
            
            % Enter the sentivity
            switch k
                case 1
                    A = Asen(:,j);
                case 2
                    B = Bsen(:,j);
                case 3
                    D = Dsen(:,j);
            end
            
            % Calculate stiffness for each element
            
            [k_el_] = FEM_elstiff(A,B,D,eltyp,crds,conn,numind,nqpts,iele);
            
            de_dxi(k,j) = -0.5*d_el'*k_el_*d_el;
            
        end
    end
    
    if typOPT == 1
        
        switch typLAM
        
            case 'ORT'
        
                 grad_obj = [grad_obj; de_dxi(1,1); de_dxi(1,3);];
                
            case 'GEN'
                
                grad_obj = [grad_obj; de_dxi(1,1:4)'];
                
        end
        
    elseif typOPT == 2
        
           switch typLAM
        
            case 'ORT'
                
                 grad_obj = [grad_obj; de_dxi(3,1); de_dxi(3,3);];
        
            case 'GEN'
                
                grad_obj = [grad_obj; de_dxi(3,1:4)'];
                
           end
        
        %     elseif typOPT == 3
        %         grad_obj = [grad_obj; de_dxi(1,1:4)'; de_dxi(3,1:4)'];
    end
    
end

return
