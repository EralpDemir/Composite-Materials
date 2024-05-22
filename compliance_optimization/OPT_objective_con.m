function [obj, grad_obj] = OPT_objective_con(x,numel, nnpe, dof, numind, nqpts, dofst, ndir, conn,crds, ...
    All_xi_A, All_xi_B, All_xi_D, t, U, dofBCs, valBCs, numBCs, f, K, ...
    qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, qL_u, Fipedlu_ElNod, qL_v, ...
    Fipedlv_ElNod, eltyp, Asen, Bsen, Dsen, typOPT, typLAM)

% Assign Stiffness parameters
B(1:6) = zeros(6,1);

% Set up the optimization variable
for iele=1:1:numel
    
    switch typLAM
    
        case 'ORT'
            
            if typOPT == 1

                All_xi_A(1,iele) = x(1);
                All_xi_A(3,iele) = x(2);

            elseif typOPT == 2

                All_xi_D(1,iele) = x(1);
                All_xi_D(3,iele) = x(2);

            end
            
            
            
        case 'GEN'
            
            
            if typOPT == 1

                All_xi_A(:,iele) = x;

            elseif typOPT == 2

                All_xi_D(:,iele) = x;

            end
            
            
    end
end

% Run FEM analysis for the new values of Lamination Parameters
[d, obj] = FEM_lam(numel, nnpe, dof, numind, nqpts, dofst, ndir, conn,...
    crds, All_xi_A, All_xi_B, All_xi_D, t, U, dofBCs, valBCs, numBCs, ...
    f, K, qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, qL_u, Fipedlu_ElNod, qL_v, ...
    Fipedlv_ElNod, eltyp);

grad_obj = [];







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

        
        fe=zeros(dof*nnpe,1);
        K_=K;
        
        % Assign stiffness for each element
        for iele=1:1:numel
            
            % Calculate stiffness for each element
            [ke] = FEM_elstiff(A,B,D,eltyp,crds,conn,numind,nqpts,iele);
            
            % Assemble elemental stiffness to global stiffness matrix
            [K_] = FEM_assemble(iele, nnpe, conn, numind, dofst, ndir, ke, fe, K_, f, dof, eltyp);
            
        end
        
        

        
                   
        de_dxi(k,j) = -0.5*d'*K_*d;

    end
end



if typOPT == 1
    
    switch typLAM
        
        case 'ORT'
             
            grad_obj = [    de_dxi(1,1); 
                            de_dxi(1,3)] ;           
            
        case 'GEN'
    
            grad_obj = de_dxi(1,1:4)';
            
    end

elseif typOPT == 2
    
     switch typLAM
        
        case 'ORT'
            
            grad_obj = [    de_dxi(3,1); 
                            de_dxi(3,3)] ;
            
        case 'GEN'
    
            grad_obj = de_dxi(3,1:4)';
            
     end

end
    


return
