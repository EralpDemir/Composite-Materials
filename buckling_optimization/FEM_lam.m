
function [R_m_bottom, lambda, a, di, mK_m, nKg] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
    conn, crds, All_xi_A, All_xi_D, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
    numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m, ...
    dofBCs_m_bottom, dofBCs_m_top)





% GLOBAL ROUTINE FOR FEM ANALYSIS OF ISOTROPIC COMPOSITE LAMINATES

%______________ Stiffness and Force vector calculations _________________


% BENDING STIFFNESSS CALCULATION
% Calculate element stiffnesses and assemble
for iele = 1:1:numel
    
    % Lamination parameter A of the element
    xi_A = All_xi_A(:,iele);   
    
    % Lamination parameter D of the element
    xi_D = All_xi_D(:,iele);
    
    % Material stiffness
    [A, D] = MAT_stiffnesslp(xi_A, xi_D, t, U);
    
    
    
    % Calculate stiffness for each element
    [k_b] = FEM_elstiff(A,D,eltyp_b,crds,conn,iele);
       
%      % Out-of-plane distributed load along z-direction over the surface of an element
%     [feS_w] = FEM_fsdlw(qS_w,Fsdlw_El,crds,conn,eltyp_b,nnpe_b,dof_b,iele);
%     
%     % Out-of-plane distributed load along z-direction over the edge of an element
%     [feL_w] = FEM_fedlw(qL_w,nnpe,Fedlw_ElNod,conn,crds,eltyp_b,dof_b,iele);
% 
%     % Total force vector
%     fe =feS_w + feL_w;
%     
    % Assemble elemental stiffness to global stiffness matrix
    [K_b] = FEM_assemble(iele,nnpe,conn,k_b,f_b,K_b,F_b,dof_b);

end


% Store stiffness matrix and force vector before modification
K_bm = K_b;
%F_bm = F_b;

% Apply Boundary conditions
%[Kbm, ~] = FEM_applybcs(Kbm, Fbm, dofBCs_b, valdofs_b, numdofs_b);
[K_bm] = FEM_modifystiffness(K_bm, dofBCs_b, numdofs_b);




% MEMBRANE STIFFNESSS CALCULATION
% Calculate element stiffnesses and assemble
for iele = 1:1:numel
    
    % Lamination parameter A of the element
    xi_A = All_xi_A(:,iele);   
    
    % Lamination parameter D of the element
    xi_D = All_xi_D(:,iele);
    
    % Material stiffness
    [A, D] = MAT_stiffnesslp(xi_A, xi_D, t, U);
    
    
    
    % Calculate stiffness for each element
    [k_m] = FEM_elstiff(A,D,eltyp_m,crds,conn,iele);
       
  
    % In-plane distributed load along x-direction over the edge of an element
    [feL_u] = FEM_fipedlu(qL_u,nnpe,Fipedlu_ElNod,crds,conn,eltyp_m,dof_m,iele);

    % In-plane distributed load along y-direction over the edge of an element
    [feL_v] = FEM_fipedlv(qL_v,nnpe,Fipedlv_ElNod,crds,conn,eltyp_m,dof_m,iele);


%      % Out-of-plane distributed load along z-direction over the surface of an element
%     [feS_w] = FEM_fsdlw(qS_w,Fsdlw_El,crds,conn,eltyp_b,nnpe_b,dof_b,iele);
%     
%     % Out-of-plane distributed load along z-direction over the edge of an element
%     [feL_w] = FEM_fedlw(qL_w,nnpe,Fedlw_ElNod,conn,crds,eltyp_b,dof_b,iele);
% 
    % Total force vector
    f_m =feL_u + feL_v;
%     
    % Assemble elemental stiffness to global stiffness matrix
    [K_m,F_m] = FEM_assemble(iele,nnpe,conn,k_m,f_m,K_m,F_m,dof_m);

end




% Store stiffness matrix and force vector before modification
mK_m=K_m;

% % This is needed to not to modify K_m (keep it as zero!)
% Km = K_m;


mF_m = F_m;

% Apply Boundary conditions
[mK_m, mF_m] = FEM_applybcs(mK_m, mF_m, dofBCs_m, valdofs_m, numdofs_m);



% IN-PLANE
% Nodal displacements
di = mK_m\mF_m;

% Reaction forces
F_m = K_m * di;
R_m_bottom = sum(F_m(dofBCs_m_bottom));
R_m_top = sum(F_m(dofBCs_m_top));



% GEOMETRIC STIFFNESSS CALCULATION
% Geometric stiffness matrix
% Calculate geometric stiffnesses and assemble


for iele = 1:1:numel

    % Lamination parameter A of the element
    xi_A = All_xi_A(:,iele);   
    
    % Lamination parameter D of the element
    xi_D = All_xi_D(:,iele);
    
    % Material stiffness
    [A] = MAT_stiffnesslp(xi_A, xi_D, t, U);
    
    % Elemental geometric stiffness matrix
    [k_g] = FEM_geostiffness(di,A,eltyp_b,crds,conn,dof_m,iele);
    
     % Assemble elemental stiffness to global stiffness matrix
    [K_g] = FEM_assemble(iele,nnpe,conn,k_g,f_g,K_g,F_g,dof_b);

end

% % Modify stiffness matrix for the given BCs
% K_g = K_g * t;

% Modify stiffness matrices to eliminate in-plane dofs
K_gm=K_g;

%[K_gm, ~] = FEM_applybcs(K_gm, F_g, dofBCs_b, valdofs_b, numdofs_b);
[K_gm] = FEM_modifystiffness(K_gm, dofBCs_b, numdofs_b);






% [Km] = FEM_modifystiffness(K, dofBCs, numdofs, numind, dof);
% [Mm] = FEM_modifystiffness(M, dofBCs, numdofs, numind, dof);

% evec=0;evals=0;
[evecs, evals] = eig(K_bm,K_gm);



 



% % Sort the eigenvalues
% if sign(min(real(diag(evals))))<0
%     [eval,ind] = sort(diag(real(evals)),'ascend');
% else
%     [eval,ind] = sort(diag(real(evals)),'descend');
% end


[eval,ind] = sort(diag(real(evals)),'descend');

eval=real(eval);






i=1;

while eval(i)==Inf || eval(i)>0
    i=i+1;
end




% % Real number (not complex)
% while ~isreal(eval(i))
%     i=i+1;
% end


lambda=eval(i);
a=evecs(:,ind(i));


a(dofBCs_b)=0;



% % Plot the mode shape
% lambda
% FEM_output_ms(a,dof_b,conn,crds,numel,nnpe,eltyp_b,1)


% % Select the minimum positive eigenvalue and corresponding eigenvector
% % The value of lambda shall be different than "1"
% lambda=1;
% i=1;
% while round(lambda) ==1
%    lambda =eval(i);
%    a =  evecs(:,ind(i));
%    i=i+1;
% end


% lambda=lambda * t;

% lambda =1;
% i=1;
% while lambda>=0
%     lambda=eval(end-i);
%     a=evecs(:,ind(end-i));
%     i = i + 1;
% end







nKg = a'*K_gm*a;



% a = a/sqrt(abs(nKg));
% 
% 
% nKg = a'*K_g*a;



% % Plot the mode shape
% lambda
% FEM_output_ms(a,dof_b,conn,crds,numel,nnpe,eltyp_b,1)





% % Plot the mode shapes
% size(eval,1)
% for i=1:size(eval,1)
%     a0=evecs(:,ind(i));
%     disp(['mode: ', num2str(i), '  frequency: ', num2str(eval(i))])
%     FEM_output_ms(a0,dof_b,conn,crds,numel,nnpe,eltyp_b,10)
% end

 
 %FEM_output(d,dof,dofst,conn,crds,numel,nnpe,eltyp,scale)

 

%___________________________ Post-Processing ___________________________

% % Reaction forces
% r = K*d-f;
% 
% % Strain Energy
% e = 1/2*d'*K*d;

% Plot the result using the element shape functions
%FEM_output(d,dof,dofst,conn,crds,numel,nnpe,eltyp,scale)

% % Plot the stress distribution and find the location, value, and related
% % elements of maximum stresses
% [max_stress_locs, max_stress_vals, max_stress_els] = ...
%     FEM_output_stressdist(nPly,numel,conn,crds,nnpe,dof,dofst,numind,...
%     d,Q_bar,angPly,zPly,eltyp)

% % Plot the mesh with the maximum stresses shown on it
% FEM_output_mesh(max_stress_locs,d,dofst,conn,...
%     crds,numel,nnpe,eltyp,scale)
