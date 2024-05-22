function [obj, grad_obj] = OPT_objective_var(x, numel, nnpe, dof, ...
    conn, crds, rho, t, U, dofBCs, valdofs, numdofs, f, F, K, M, ...
    eltyp, Dsen, typLAM)



All_xi_A = zeros(4,numel);
All_xi_D = zeros(4,numel);


switch typLAM
    
    case 'ORT'

        
        All_xi_D(1,1:numel) = x(1:2:end);
        All_xi_D(3,1:numel) = x(2:2:end);
        
        
        
    case 'GEN'

        for i=1:1:numel
            All_xi_D(1,i) = x(4*i-3,1);
            All_xi_D(2,i) = x(4*i-2,1);
            All_xi_D(3,i) = x(4*i-1,1);
            All_xi_D(4,i) = x(4*i,1);
        end

end






% Run FEM analysis for the new values of Lamination Parameters
[lambda, a, nM] = FEM_lam(numel, nnpe, dof, conn, crds, rho, ...
    All_xi_A, All_xi_D, t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp);


% if lambda<0
%     mult = 1;
% else
%     mult = -1;
% end
% 
% 
% obj =  mult * lambda;

obj =  -lambda ;














% Calculate Bending Stiffness
dkb_dxi = zeros(numel,2);





for iele = 1:1:numel
    
    % Elemental displacements
    a_i = [];
    for j = 1:1:nnpe
        a_i = [a_i; a(dof*conn(iele,j)-dof+1:1:dof*conn(iele,j))];
    end
    
    % Compute the sensitivities:
    % Bending stiffness
    % With respect only the lamination parameter D
    for k = 1:1:4
        
    
            % Assign Stiffness parameters

            A = zeros(6,1);    
            D = zeros(6,1);
            
            % Enter the sentivity
            D = Dsen(:,k);
           
            
            % Calculate stiffness for each element
                        
            
            [k_bsen] = FEM_elstiff(A,D,eltyp,crds,conn,iele);
            
            
            % Local sensitivity of bending stiffness
            dkb_dxi(iele,k) = a_i' * k_bsen * a_i;
            


    end
    
    
    

    
    
    
    
    
    
end
% End of the elemental loop










    
grad_obj = [];

for iele = 1:1:numel    


     switch typLAM
        
        case 'ORT'
    
            grad_obj = [   grad_obj; 
                dkb_dxi(iele,1); 
                dkb_dxi(iele,3); ];
            
            
            
        case 'GEN'
           
            grad_obj = [   grad_obj; 
                    dkb_dxi(iele,1); 
                    dkb_dxi(iele,2); 
                    dkb_dxi(iele,3); 
                    dkb_dxi(iele,4);];

                
     end
                
end



grad_obj =  -grad_obj/nM;




return
