function [d, e, r] = FEM_lam(numel, nnpe, dof, numind, nqpts, dofst, ...
    ndir, conn, crds, All_xi_A, All_xi_B, All_xi_D, t, U, dofBCs, ...
    valBCs, numBCs, f, K, qS_w, Fsdlw_El, qL_w, Fedlw_ElNod, ...
    qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp)

%FEM_lam     calculates the global routine for FEM analysis of
%            composite laminates.


%_________________ Stiffness and Force vector calculations _______________
% Calculate element stiffnesses and assemble
for iele = 1:1:numel
    
    % Store all of the lamination parameter PER ELEMENT
    xi_A = All_xi_A(1:4,iele);
    xi_B = All_xi_B(1:4,iele);
    xi_D = All_xi_D(1:4,iele);
    
    % A matrix components - lamina
    A = t*U*[1; xi_A] ;
    
    % B matrix components - lamina
    B = t^2/4*U*[0; xi_B];
    
    % D matrix components - lamina
    D = t^3/12*U*[1; xi_D] ;
    
    % Calculate stiffness for each element
    [ke] = FEM_elstiff(A,B,D,eltyp,crds,conn,numind,nqpts,iele);
    
    % In-plane distributed load along x-direction over the edge of an element
    [feL_u] = FEM_fipedlu(qL_u,nnpe,Fipedlu_ElNod,crds,conn,eltyp,dof,iele,numind);
    
    % In-plane distributed load along y-direction over the edge of an element
    [feL_v] = FEM_fipedlv(qL_v,nnpe,Fipedlv_ElNod,crds,conn,eltyp,dof,iele,numind);
    
    % Out-of-plane distributed load along z-direction over the surface of an element
    [feS_w] = FEM_fsdlw(qS_w,Fsdlw_El,crds,conn,eltyp,nnpe,dof,iele,numind,nqpts);
    
    % Out-of-plane distributed load along z-direction over the edge of an element
    [feL_w] = FEM_fedlw(qL_w,nnpe,Fedlw_ElNod,conn,crds,eltyp,dof,iele,numind);
    
    % Total force vector
    fe = feL_u + feL_v + feS_w + feL_w;
    
    % Assemble elemental stiffness to global stiffness matrix
    [K,f,ndir] = FEM_assemble(iele, nnpe, conn, numind, dofst, ndir, ke, fe, K, f, dof, eltyp);
    
end

%_______________________ finding the displacement _______________________

% Store stiffness matrix and force vector before modification

Km = K;
fm = f;

% Apply Boundary conditions
[Km, fm] = FEM_applybcs(Km, fm, dofBCs, valBCs, numBCs);

% Nodal displcements
d = Km\fm;

% Reaction forces
r = K*d - f;

% Strain Energy
e = 1/2*d'*K*d;

% % Plot the result using the element shape functions
% FEM_output(d,dof,dofst,conn,crds,numel,nnpe,eltyp,scale)

