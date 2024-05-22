function [obj, grad_obj] = OPT_objective_var(x, numel, nnpe, dof_b, dof_m, ...
    conn, crds, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m,f_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, numind_m, ...
    numind_b, Asen, Dsen, BigNmat, Bigfvec, LSC_nnpe, ...
    LSC_conn, LSC_totnod, N0,dofBCs_m_bottom, dofBCs_m_top)


All_xi_A = zeros(4,numel);
All_xi_D = zeros(4,numel);









% Smooth the lamination parameters
[x_nodes] = LSC_LP(BigNmat, Bigfvec, numel, LSC_nnpe, ...
    LSC_conn, LSC_totnod, x);

% Calculate the angles at the element centers
[x_c] = LSC_LP0(numel, LSC_nnpe, LSC_conn, ...
    x_nodes, N0);



for iele=1:1:numel
    All_xi_A(1,iele) = x_c(iele,1);
    All_xi_A(3,iele) = x_c(iele,2);
    All_xi_D(1,iele) = x_c(iele,3);
    All_xi_D(3,iele) = x_c(iele,4);
end


% Run FEM analysis for the new values of Lamination Parameters
[R_m_bottom,lambda, a, di, mK_m, nKg] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
    conn, crds, All_xi_A, All_xi_D, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m,...
    dofBCs_m_bottom, dofBCs_m_top);



if lambda < 0
    
    mult = 1;
    
else
    
    mult = -1;
    
end
% 
% 

obj = mult *lambda ;

% Initialization of force vector and stiffness matrices


grad_obj = [];



% Compute the sensitivities:



dkb_dxi = zeros(numel,2);



% Bending stiffness

% Assign stiffness for each element
for iele=1:1:numel

    
    % Elemental displacements
    ai = zeros(nnpe*dof_b,1);
    for j = 1:1:nnpe
        ai((j-1)*dof_b+1:1:j*dof_b )=  a(dof_b*conn(iele,j)-dof_b+1:1:dof_b*conn(iele,j));
    end
    
    

    % With respect only the lamination parameter D
    for k = 1:1:2

        % Assign Stiffness parameters
        A = zeros(6,1);    
        D = zeros(6,1);

        % Enter the sentivity
        switch k

            case 1
                D = Dsen(:,1);
            case 2
                D = Dsen(:,3);
        end



    
    

        
%         % Initialize the sensitivty of bending stiffness
%         K_bsen =zeros(numind_b,numind_b);
        

        % Calculate stiffness for each element
        [k_bsen] = FEM_elstiff(A,D,eltyp_b,crds,conn,iele);

%         % Assemble elemental stiffness to global stiffness matrix
%         [K_bsen] = FEM_assemble(iele,nnpe,conn,k_bsen,f_b,K_bsen,F_b,dof_b);
% 
%         
%         
%         % Enter the sentivity
%         dkb_dxi(iele,k) = a' * K_bsen * a;
        
      
        
        dkb_dxi(iele,k) = ai' * k_bsen * ai;
        
        
        
    end
    
    
%     [K_bsen] = FEM_modifystiffness0(K_bsen, dofBCs_b, numdofs_b); 
    
    
        


   
                   



end





    
% Calculation of Geometric stiffness
% 1st derivative
% With respect only the lamination parameter A

% Initialize

ddi_dxiA1 = zeros(numind_m,numel);
ddi_dxiA3 = zeros(numind_m,numel);
dkg1_dxi = zeros(numel,2);



% Assign stiffness for each element
for iele=1:1:numel

    
    
    % Elemental displacements
    ai = zeros(nnpe*dof_b,1);
    for j = 1:1:nnpe
        ai((j-1)*dof_b+1:1:j*dof_b )=  a(dof_b*conn(iele,j)-dof_b+1:1:dof_b*conn(iele,j));
    end
    
    
    
    
    for k = 1:1:2



        % Assign Stiffness parameters
        A = zeros(6,1);
        D = zeros(6,1);


        % Enter the sentivity
        switch k
            case 1
                A = Asen(:,1);
            case 2
                A = Asen(:,3);
        end


    
    


        % Calculate geometric stiffness for each element
        [k_gsen1] = FEM_geostiffness(di,A,eltyp_b,crds,conn,dof_m,iele);
        
         % Membrane stiffness
        [k_msen] = FEM_elstiff(A,D,eltyp_m,crds,conn,iele);



        
        dKm_dxiA1 = zeros(numind_m,numind_m);
        dKm_dxiA3 = zeros(numind_m,numind_m);
        
        
        
%          dkg1_dxi(iele,k) = -lambda * a' * K_gsen1 * a;
         
         dkg1_dxi(iele,k) = -lambda * ai' * k_gsen1 * ai;
        
        
        % Assemble membrane stiffness
        switch k

            case 1
                
                [dKm_dxiA1] = FEM_assemble(iele,nnpe,conn,k_msen,f_m,dKm_dxiA1,F_m,dof_m);

                
                % Derivatives of in-plane displacements
                ddi_dxiA1(1:numind_m,iele)= -inv(mK_m) * dKm_dxiA1 * di;
                
                
                
                
                
            case 2
        
                [dKm_dxiA3] = FEM_assemble(iele,nnpe,conn,k_msen,f_m,dKm_dxiA3,F_m,dof_m);
                
               
                
                % Derivatives of in-plane displacements
                ddi_dxiA3(1:numind_m,iele)= -inv(mK_m) * dKm_dxiA3 * di;
        end
                

        
        
        
        
                
    end
    
    
    
    
    
    
    

    


end










% The non-local term
dkg2_dxi = zeros(numel,2);


for k = 1:1:2
    
    
        


        
        switch k

            case 1
                
                
                for iel=1:numel
                
                
                
%                 A = Asen(:,1);

                    ddi(1:numind_m,1) = ddi_dxiA1(1:numind_m,iel);
                    
                    % Set the derivatives of prescribed displacements to ZERO!
                    [ddi] = OPT_modifydisp(ddi, dofBCs_m, numdofs_m);
                
                
                    dKg2_dxiA1 = zeros(numind_b,numind_b);
                    
                    for iele=1:numel
                
                        
                    
                    
                    
                        % Lamination parameter A of the element
                        xi_A = All_xi_A(:,iele);   

                        % Lamination parameter D of the element
                        xi_D = All_xi_D(:,iele);

                        % Material stiffness
                        [A] = MAT_stiffnesslp(xi_A, xi_D, t, U);






                        % Calculate stiffness for each element
                        [k_gsen2] = FEM_geostiffness(ddi,A,eltyp_b,crds,conn,dof_m,iele);


                        % Assemble
                        [dKg2_dxiA1] = FEM_assemble(iele,nnpe,conn,k_gsen2,f_g,dKg2_dxiA1,F_g,dof_b);

                     
                   
                     
                    end
                 
                    
                    
                    
                    dkg2_dxi(iel, k) = -lambda * a' * dKg2_dxiA1 * a;
                 
               

                
                
                end
                
                
                
            case 2

%                 A = Asen(:,3);


                for iel=1:numel
                
                    ddi(1:numind_m,1) = ddi_dxiA3(1:numind_m,iel);
                
                    
                    
                    % Set the derivatives of prescribed displacements to ZERO!
                    [ddi] = OPT_modifydisp(ddi, dofBCs_m, numdofs_m);

                    
                    
                    dKg2_dxiA3 = zeros(numind_b,numind_b);
                 
                    % Assemble 
                    for iele=1:1:numel

                       

                        % Lamination parameter A of the element
                        xi_A = All_xi_A(:,iele);   

                        % Lamination parameter D of the element
                        xi_D = All_xi_D(:,iele);

                        % Material stiffness
                        [A] = MAT_stiffnesslp(xi_A, xi_D, t, U);                    


                        % Calculate stiffness for each element
                        [k_gsen2] = FEM_geostiffness(ddi,A,eltyp_b,crds,conn,dof_m,iele);


                        % Assemble
                        [dKg2_dxiA3] = FEM_assemble(iele,nnpe,conn,k_gsen2,f_g,dKg2_dxiA3,F_g,dof_b);




                    end
                 
                    
                    
                    dkg2_dxi(iel, k) = -lambda* a' * dKg2_dxiA3 * a;
                  
                    
       
                 
                end


        end
        
        

        
        

        
        
end


% Derivative of bending stiffness
dkg_dxi =   dkg1_dxi + dkg2_dxi;

grad_obj = [];

for iele = 1:1:numel    


           
     grad_obj = [   grad_obj;
                    dkg_dxi(iele,1); 
                    dkg_dxi(iele,2); 
                    dkb_dxi(iele,1); 
                    dkb_dxi(iele,2) ];

                
end


grad_obj =  grad_obj/nKg * mult ;





return
