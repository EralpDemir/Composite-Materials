


% Initialize the size of the matrices
BigMmat = zeros(totnod,totnod);






% Total area of the domain
Area = 0;

% Matrices for interpolation
for iele = 1:1:numel
    
    % Nodes from connectivity
    nodes = conn(iele,:);
    
    % Element coordinates
    elcrds=crds(nodes,1:2);
    
%     % Angles at the nodes of the element (radians)
%     angs(1:nnpe,1) = angles(nodes);
    
%     % Angle - constant throughout the element
%     ang = mean(angs);
    
    % Shape function derivatives and determinant
    [dndx,dndy,detj] = jacobian(nnpe, nqpt, dN_qpt, elcrds);
    
%     % Least Square mapping and pseudo force vector
%     [Mmat, fvec, area] = LSC_mmat(LSC_nnpe, LSC_nqpt, LSC_N_qpt, LSC_wght, detj);

    % Coefficient matrix and pseudo-force vector
%     [Mmat, Mvec, area] = mmat(nnpe, nqpt, N_qpt, wght, detj, angs);
    [Mmat, Mvec, area] = mmat(nnpe, nqpt, N_qpt, wght, detj, angs);
    
    
    % Assemble the local elemental constraints to the corresponding global
    % degrees of freedom
    [BigMmat] = assemblecons(iele,nnpe,conn,Mmat,BigMmat);
    
    % Assemble the pseudo-force vector
    [BigMvec] = assembleforc(iele,nnpe,conn,Nvec,BigMvec);
    

    % Overall area of the domain
    Area = Area + area;
    
end

% Interpolation functions at the center of the element
N0 = shape_fun_c(eltyp);