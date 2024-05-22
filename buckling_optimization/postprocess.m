
% % Post-processed variables
% All_angPly_post=zeros(numel,nPlys);
% 
% % Find the average of the angles
% for i=1:nPlys
%     av=mean(All_angPly(:,i));
% 
% 
%     if av>45 
%         
%         ang_range = 0:div:180-div;
%         
% 
%     elseif av<-45
% 
%         ang_range = 0:div:180-div;
%     else
% 
% 
%         ang_range = -90:div:90-div;
%     end
%     
%     disp(['average angle of layer ', num2str(i), '  :', num2str(av)])
% 
%     % Construct the Lamination Parameter Space
%     [xi_As, xi_Bs, xi_Ds, noAngPoss, angPlys] = ...
%         LP_space(ang_range, t, Ltype, nPlys, zPly);
% 
% 
%     
%     
%     
%     All_angPly_p = zeros(numel,nPlys);
%     All_xi_A_lib_p = zeros(4,numel);
%     All_xi_D_lib_p = zeros(4,numel);
% 
%     % Find the nearest lamination parameter in the space
%     for iele = 1:1:numel
% 
%         xi = [All_xi_A(1,iele); All_xi_A(3,iele); All_xi_D(1,iele); All_xi_D(3,iele);];
% 
%         [ang, xi_A, xi_D] = LP_find(xi_As, xi_Ds, angPlys, noAngPoss, xi);
% 
%         All_angPly_p(iele,1:nPlys) = ang;
% 
%         All_xi_A_lib_p(1:4,iele) = xi_A;
% 
%         All_xi_D_lib_p(1:4,iele) = xi_D;
% 
%     end    
%     
%     
%     
%     All_angPly_post(:,i) = All_angPly_p(:,i);
%     
%     
%     
% 
% end


LSC_inputs


LSC_initialize




% Smooth the angles
[angPlyNodes] = LSC_angdist(BigNmat, Bigfvec, numel, LSC_nnpe, ...
    LSC_conn, LSC_totnod, nPlys, All_angPly);

% Calculate the angles at the element centers
[All_angPly_c] = LSC_angel0(numel, LSC_nnpe, LSC_conn, nPlys, ...
    angPlyNodes, N0);

% Compute overall curvature
kappa = zeros(1,nPlys);
for iPly = 1:1:nPlys
    angs(:,1) = angPlyNodes(:,iPly)*pi/180;
    kappa(iPly) = sqrt(angs'*BigCmat*angs/Area)/2;
end

clear angs

% Compute curvature elementally
kappa_el = zeros(numel,nPlys);
for iele = 1:1:numel
    for iPly = 1:1:nPlys
        nodes = LSC_conn(iele,:);
        angs = angPlyNodes(nodes,iPly)*pi/180;
        Cmat(:,:) = Cmat_els(iele,:,:);
        kappa_el(iele,iPly) = sqrt(angs'*Cmat*angs/sx/sy)/2;
    end
end

clear angs nodes iele iPly



% 
% % Run FEM analysis for the final/optimized values of Lamination Parameters
% [lambda_lib, a_lib] = FEM_lam(numel, nnpe, dof_b, dof_m, ...
%     conn, crds, All_xi_A_lib_p, All_xi_D_lib_p, t, U, dofBCs_b, dofBCs_m, valdofs_b, valdofs_m, ...
%     numdofs_b, numdofs_m, f_b, F_b, K_b, F_m, K_m, f_g, F_g, K_g, ...
%     qL_u, Fipedlu_ElNod, qL_v, Fipedlv_ElNod, eltyp_b, eltyp_m);

    
    