


% Initialize the size of the matrices
BigNmat = zeros(totnod,totnod);
BigCmat = zeros(totnod,totnod);
BigDmat = zeros(totnod,totnod);






Bigfvec = zeros(totnod,1);


% Total area of the domain
Area = 0;

% Matrices for interpolation
for iele = 1:1:numel
    
    % Nodes from connectivity
    nodes = conn(iele,:);
    
    % Element coordinates
    elcrds=crds(nodes,1:2);
    
    
    
    
%     angs = angles(iele)*inv(N_qpt)*ones(nnpe,1);
%     
    
    % Angles at the nodes of the element (radians)
    angs(1:nnpe,1) = angles(nodes);
        
    
    
%     % Angle - constant throughout the element
%     ang = mean(angs);
    
    % Shape function derivatives and determinant
    [dndx,dndy,detj] = jacobian(nnpe, nqpt, dN_qpt, elcrds);
    
%     % Least Square mapping and pseudo force vector
%     [Mmat, fvec, area] = LSC_mmat(LSC_nnpe, LSC_nqpt, LSC_N_qpt, LSC_wght, detj);

    % Coefficient matrix and pseudo-force vector
    [Nmat, Cmat, Dmat, Nvec, area] = nmat(nnpe, nqpt, wght, N_qpt, dndx, dndy, detj, angs);
    
    % Assemble the local elemental constraints to the corresponding global
    % degrees of freedom
    [BigNmat] = assemblecons(iele,nnpe,conn,Nmat,BigNmat);
    [BigCmat] = assemblecons(iele,nnpe,conn,Cmat,BigCmat);
    [BigDmat] = assemblecons(iele,nnpe,conn,Dmat,BigDmat);
    
    % Assemble the pseudo-force vector
    [Bigfvec] = assembleforc(iele,nnpe,conn,Nvec,Bigfvec);

    % Overall area of the domain
    Area = Area + area;
    
end

% Interpolation functions at the center of the element
N0 = shape_fun_c(eltyp);