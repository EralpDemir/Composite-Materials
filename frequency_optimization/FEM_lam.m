
function [lambda, a, nM] = FEM_lam(numel, nnpe, dof, conn, crds, rho, ...
    All_xi_A, All_xi_D, t, U, dofBCs, valdofs, numdofs, f, F, K, M, eltyp)






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
    [k] = FEM_elstiff(A,D,eltyp,crds,conn,iele);
       
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
    [K] = FEM_assemble(iele,nnpe,conn,k,f,K,F,dof);

end


% Store stiffness matrix and force vector before modification
K_m = K;
%F_bm = F_b;

% Apply Boundary conditions
[K_m, ~] = FEM_applybcs(K_m, F, dofBCs, valdofs, numdofs);
%[K_m] = FEM_modifystiffness(K_m, dofBCs, numdofs);



% MASS MATRIX CALCULATION

for iele = 1:1:numel

    
    % Elemental geometric stiffness matrix
    [m] = FEM_massmatrix(rho,t,eltyp,crds,conn,iele);
    
     % Assemble elemental stiffness to global stiffness matrix
    [M] = FEM_assemble(iele,nnpe,conn,m,f,M,F,dof);

end



% Modify stiffness matrices to eliminate in-plane dofs
M_m=M;

%[M_m, ~] = FEM_applybcs(M_m, F, dofBCs, valdofs, numdofs);
[M_m] = FEM_modifystiffness(M_m, dofBCs, numdofs);






% [Km] = FEM_modifystiffness(K, dofBCs, numdofs, numind, dof);
% [Mm] = FEM_modifystiffness(M, dofBCs, numdofs, numind, dof);

% evec=0;evals=0;
[evecs, evals] = eig(K_m,M_m);



     




% Sort the eigenvalues

[eval,ind] = sort(diag(real(evals)),'ascend');



lambda=eval(1);
a=evecs(:,ind(1));



% 
% 
% 
% % Sort the eigenvalues
% [eval,ind] = sort(diag(evals),'ascend');
% 
% 
% % lambda=eval(1);
% % a=evecs(:,ind(1));
% 
% 
% 
% % Select the minimum positive eigenvalue and corresponding eigenvector
% % The value of lambda shall be different than "1"
% lambda=1;
% i=1;
% while round(lambda) == 1 
%    lambda =eval(i);
%    a =  evecs(:,ind(i));
%    i=i+1;
% end






nM = a'*M*a;


% % Normalize eigenvectors s.t. a' * M * a = 1
% a = a/sqrt(abs(nM));
% 
% 
% 
% 
% 
% % check normalization 
% nM = a'*M*a;






% % Plot the mode shape
% lambda
% FEM_output_ms(a,dof,conn,crds,numel,nnpe,eltyp,1)





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
return
