function [obj, grad_obj] = OPT_objective_con(x, numel, nnpe, dof_b, dof_m, ...
    conn, crds, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m,f_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, numind_m, ...
    numind_b, Asen, Dsen, dofBCs_m_bottom, dofBCs_m_top)

% Assign Stiffness parameters
All_xi_A = zeros(4,numel);
All_xi_D = zeros(4,numel);

% Set up the optimization variable
for iele=1:1:numel
    

            


    All_xi_A(1,iele) = x(1);
    All_xi_A(3,iele) = x(2);

    All_xi_D(1,iele) = x(3);
    All_xi_D(3,iele) = x(4);

            
end


% Run FEM analysis for the new values of Lamination Parameters
[R_m_bottom, lambda, a, di, mK_m, nKg] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
    conn, crds, All_xi_A, All_xi_D, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m,...
    dofBCs_m_bottom, dofBCs_m_top);



if lambda<0%* R_m_bottom < 0
    
    mult = 1;
    
else
    
    mult = -1;
    
end




obj = mult*lambda;



% Initialization of force vector and stiffness matrices


grad_obj = [];





% Compute the sensitivities:

% Bending stiffness
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



    % Initialize the sensitivty of bending stiffness
    K_bsen =zeros(numind_b,numind_b);
    
    % Assign stiffness for each element
    for iele=1:1:numel

        % Calculate stiffness for each element
        [k_bsen] = FEM_elstiff(A,D,eltyp_b,crds,conn,iele);

        % Assemble elemental stiffness to global stiffness matrix
        [K_bsen] = FEM_assemble(iele,nnpe,conn,k_bsen,f_b,K_bsen,F_b,dof_b);

    end
    

        
    % Enter the sentivity
    switch k

        case 1

            dKb_dxiD1 = a' * K_bsen * a;
            
        case 2

            dKb_dxiD3 = a' * K_bsen * a;

    end

   
                   



end





    
% Calculation of Geometric stiffness
% 1st derivative
% With respect only the lamination parameter A

% Initialize
dKm_dxiA1 = zeros(numind_m,numind_m);
dKm_dxiA3 = zeros(numind_m,numind_m);



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


    K_gsen1=zeros(numind_b,numind_b);
    
    % Assign stiffness for each element
    for iele=1:1:numel

        % Calculate geometric stiffness for each element
        [k_gsen1] = FEM_geostiffness(di,A,eltyp_b,crds,conn,dof_m,iele);
        
         % Membrane stiffness
        [k_m] = FEM_elstiff(A,D,eltyp_m,crds,conn,iele);


        % Assemble elemental stiffness to global stiffness matrix
        [K_gsen1] = FEM_assemble(iele,nnpe,conn,k_gsen1,f_g,K_gsen1,F_g,dof_b);


        % Assemble membrane stiffness
        switch k

            case 1
                
                [dKm_dxiA1] = FEM_assemble(iele,nnpe,conn,k_m,f_m,dKm_dxiA1,F_m,dof_m);
        
                
            case 2
        
                [dKm_dxiA3] = FEM_assemble(iele,nnpe,conn,k_m,f_m,dKm_dxiA3,F_m,dof_m);
                
        end
                
                
    end

    % Enter the sentivity
    switch k

        case 1

            dKg1_dxiA1 = -lambda * a' * K_gsen1 * a;


        case 2

            dKg1_dxiA3 = -lambda * a' * K_gsen1 * a;

    end

    
    


end



% Derivatives of in-plane displacements
ddi_dxiA1= -inv(mK_m) * dKm_dxiA1 * di;

ddi_dxiA3= -inv(mK_m) * dKm_dxiA3 * di;




% In-plane stiffness

% Lamination parameter A of the element
xi_A = All_xi_A(:,1);   

% Lamination parameter D of the element
xi_D = All_xi_D(:,1);

% Material stiffness
[A] = MAT_stiffnesslp(xi_A, xi_D, t, U);



for k = 1:1:2
    
    
        


        
        switch k

            case 1
                
%                 A = Asen(:,1);

                ddi = ddi_dxiA1;
                
                % Set the derivatives of prescribed displacements to ZERO!
                [ddi] = OPT_modifydisp(ddi, dofBCs_m, numdofs_m);
                
                
            case 2

%                 A = Asen(:,3);
                
                ddi = ddi_dxiA3;
                
                % Set the derivatives of prescribed displacements to ZERO!
                [ddi] = OPT_modifydisp(ddi, dofBCs_m, numdofs_m);

        end
        
        
        K_gsen2=zeros(numind_b,numind_b);
        
        
        
        
        
         % Assign stiffness for each element
        for iele=1:1:numel
            
            % Calculate stiffness for each element
            [k_gsen2] = FEM_geostiffness(ddi,A,eltyp_b,crds,conn,dof_m,iele);
        
        
            % Assemble elemental stiffness to global stiffness matrix
            [K_gsen2] = FEM_assemble(iele,nnpe,conn,k_gsen2,f_g,K_gsen2,F_g,dof_b);
            
        end
        
        
     
        
        switch k

            case 1
                

               dKg2_dxiA1 = -lambda* a'*K_gsen2*a;
                
                
            case 2

                
               dKg2_dxiA3 = -lambda* a'*K_gsen2*a;

        end
        
        
        
end


% Derivative of bending stiffness
dobj_dxiD1 =  dKb_dxiD1 ;

dobj_dxiD3 =  dKb_dxiD3 ;


% Derivative of geometric stiffness
dobj_dxiA1 =  (dKg1_dxiA1 + dKg2_dxiA1);

dobj_dxiA3 =  (dKg1_dxiA3 + dKg2_dxiA3);


    
grad_obj =  [  dobj_dxiA1; 
               dobj_dxiA3;
               dobj_dxiD1;
               dobj_dxiD3;] * mult / nKg ;%  * R_m_bottom;           




 


return
