function [obj, grad_obj] = OPT_objective_con(x, numel, nnpe, dof, ...
    conn, crds, rho, t, U, dofBCs, valdofs, numdofs, f, F, K, M, ...
    eltyp, Dsen, typLAM)

% Assign Stiffness parameters

All_xi_A = zeros(4,numel);
All_xi_D = zeros(4,numel);

% Set up the optimization variable

 switch typLAM
    
        case 'ORT'
            
            
            for iele=1:1:numel
               All_xi_D(1,iele) = x(1);


               All_xi_D(3,iele) = x(2);

            end
            
            
        case 'GEN'    

            for iele=1:1:numel



                All_xi_D(1,iele) = x(1);
                All_xi_D(2,iele) = x(2);

                All_xi_D(3,iele) = x(3);
                All_xi_D(4,iele) = x(4);


            end
 end
 

% Run FEM analysis for the new values of Lamination Parameters
[lambda, a, nM] = FEM_lam(numel, nnpe, dof, conn, crds, rho, ...
    All_xi_A, All_xi_D, t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp);



% if lambda < 0
%     
%     mult = 1;
%     
% else
%     
%     mult = -1;
%     
% end
% 
% 
% obj = mult * lambda;

obj = - lambda ;

% Initialization of force vector and stiffness matrices



dKb_dxiD=zeros(4,1);

% Bending stiffness
% With respect only the lamination parameter D
for k = 1:1:4

    % Assign Stiffness parameters
    A = zeros(6,1);    
    D = zeros(6,1);

    % Enter the sentivity
    D = Dsen(:,k);





    % Initialize the sensitivty of bending stiffness
    K_bsen =K;
    
    % Assign stiffness for each element
    for iele=1:1:numel

        % Calculate stiffness for each element
        [k_bsen] = FEM_elstiff(A,D,eltyp,crds,conn,iele);

        % Assemble elemental stiffness to global stiffness matrix
        [K_bsen] = FEM_assemble(iele,nnpe,conn,k_bsen,f,K_bsen,F,dof);

    end
    
    
    
    
        
     % Enter the sentivity
      dKb_dxiD(k,1) = a' * K_bsen * a;
            


   
                   



end



% Derivative of bending stiffness



switch typLAM
        
        case 'ORT'
            
            grad_obj = [    dKb_dxiD(1) / nM ; 
                            dKb_dxiD(3) / nM ;] ;
            
        case 'GEN'
    
            grad_obj = - dKb_dxiD / nM ;     
            
             
      
end

    


return
