% This script initializes the LSTCC requires inputs

%_____________ Find connectivity and node coordinates ___________________
switch LSC_eltyp
    
    case 3
        
        LSC_conn = conn;
        
        LSC_crds = crds;
        
        LSC_nnpe = 3;
        
        LSC_nqpt = 1;
        
        LSC_qpt = [ 1/3 1/3];
        
        LSC_wght = 1/2;
        
    case 4
        
        LSC_conn = conn;
        
        LSC_crds = crds;
        
        LSC_nnpe = 4;
        
        LSC_qpt = [ -1/sqrt(3) -1/sqrt(3);
            1/sqrt(3) -1/sqrt(3);
            1/sqrt(3)  1/sqrt(3);
            -1/sqrt(3)  1/sqrt(3)];
        
        LSC_nqpt = 4;
        
        LSC_wght = [ 1;
            1;
            1;
            1];
        
    case 8
        
        LSC_conn = zeros(numel,8);
        
        LSC_nnpe = 8;
        
        LSC_qpt = [  -sqrt(0.6) -sqrt(0.6);
            -sqrt(0.6) 0;
            -sqrt(0.6) sqrt(0.6);
            0 -sqrt(0.6);
            0 0;
            0 sqrt(0.6);
            sqrt(0.6) -sqrt(0.6);
            sqrt(0.6) 0;
            sqrt(0.6) sqrt(0.6)];
        
        LSC_nqpt = 9;
        
        LSC_wght = [ 5/9;
            5/9;
            5/9;
            5/9;
            8/9;
            5/9;
            5/9;
            5/9;
            5/9];
        
        dum = numnod;
        
        LSC_crds=crds;
        
        for i = 1:1:numel
            
            for j = nnpe+1:1:LSC_nnpe
                
                LSC_conn(i,j-nnpe) = conn(i,j-nnpe);
                
                dum = dum+1;
                
                % Edge number
                switch j-nnpe
                    
                    case  1
                        nodeA = conn(i,1);
                        nodeB = conn(i,2);
                        x = mean([crds(nodeA,1),crds(nodeB,1)]);
                        y = mean([crds(nodeA,2),crds(nodeB,2)]);
                        
                    case 2
                        nodeA = conn(i,2);
                        nodeB = conn(i,3);
                        x = mean([crds(nodeA,1),crds(nodeB,1)]);
                        y = mean([crds(nodeA,2),crds(nodeB,2)]);
                        
                    case 3
                        nodeA = conn(i,3);
                        nodeB = conn(i,4);
                        x = mean([crds(nodeA,1),crds(nodeB,1)]);
                        y = mean([crds(nodeA,2),crds(nodeB,2)]);
                        
                    case 4
                        nodeA = conn(i,4);
                        nodeB = conn(i,1);
                        x = mean([crds(nodeA,1),crds(nodeB,1)]);
                        y = mean([crds(nodeA,2),crds(nodeB,2)]);
                        
                end
                LSC_conn(i,j) = dum;
                
                LSC_crds = [ LSC_crds;
                    x y];
                
            end
            
        end
        
        % Find coincident and indepedent nodes
        sweep = []; indep = [];
        dum = 0; dumm = 0;
        for i = 1:1:numel
            for j = nnpe+1:1:LSC_nnpe
                nodeA = LSC_conn(i,j);
                xy = LSC_crds(nodeA,:);
                entry = 0;
                for k = 1:1:numel
                    for l = nnpe+1:1:LSC_nnpe
                        nodeB = LSC_conn(k,l);
                        xy_ = LSC_crds(nodeB,:);
                        if norm(xy-xy_) < 1e-10
                            if not(i==k)
                                entry = 1;
                                dum = dum+1;
                                sweep(dum,:) = [nodeA nodeB];
                                ind(dum) = 0;
                            end
                        end
                    end
                end
                
                if entry == 0
                    dumm = dumm+1;
                    indep(dumm) = nodeA;
                end
            end
        end
        
        LSC_crds_new(1:numnod,1:2) = crds;
        LSC_conn_new = zeros(numel,LSC_nnpe);
        LSC_conn_new(1:numel,1:nnpe) = conn;
        
        % Clean up and Renumber duplicate nodes in the mesh
        no = numnod;
        
        % Add the indepedent nodes
        for i = 1:1:dumm
            no = no+1;
            nodeA = indep(i);
            LSC_crds_new(no,1:2) = LSC_crds(nodeA,1:2);
            for k = 1:1:numel
                for l = nnpe+1:1:LSC_nnpe
                    if LSC_conn(k,l) == nodeA
                        LSC_conn_new(k,l) = no;
                    end
                end
            end
        end
        
        for i = 1:1:dum
            
            nodeA = sweep(i,1);
            nodeB = sweep(i,2);
            if ind(i) == 0
                no = no+1;
                ind(i)=1;
                for j = 1:1:dum
                    if sweep(j,:) == [nodeB nodeA]
                        ind(j) = 1;
                    end
                end
                
                LSC_crds_new(no,1:2) = LSC_crds(nodeA,1:2);
                
                for k = 1:1:numel
                    for l = nnpe+1:1:LSC_nnpe
                        if LSC_conn(k,l) == nodeB || LSC_conn(k,l) == nodeA
                            LSC_conn_new(k,l) = no;
                        end
                    end
                end
                
            end
            
        end
        
        LSC_crds = LSC_crds_new;
        
        LSC_conn = LSC_conn_new;
        
        clear nodeA node B i j k l dum x y xy xy_ sweep ind
        clear dumm indep no LSC_conn_new LSC_crds_new
        
end

% Total number of nodes in the mesh
LSC_totnod = size(LSC_crds,1);

% Calculation of shape functions and their derivates at the quadrature
% points
LSC_N_qpt = zeros(LSC_nqpt,LSC_nnpe);
LSC_dN_qpt = zeros(LSC_nqpt,2,LSC_nnpe);
% LSC_dN2_qpt = zeros(LSC_nqpt,2,LSC_nnpe);

for iqpt = 1:1:LSC_nqpt
    
    k = LSC_qpt(iqpt,1);
    
    e = LSC_qpt(iqpt,2);
    
    LSC_N_qpt(iqpt,1:LSC_nnpe) = LSC_n(k,e,LSC_eltyp);
    
    LSC_dN_qpt(iqpt,1:2,1:LSC_nnpe) = LSC_nder(k,e,LSC_eltyp);
    
%     LSC_dN2_qpt(iqpt,1:2,1:LSC_nnpe) = LSC_nder2(k,e,LSC_eltyp);
    
    
end

% Initialize the size of the matrices
BigNmat = zeros(LSC_totnod,LSC_totnod);

BigCmat = zeros(LSC_totnod,LSC_totnod);

% BigDmat = zeros(LSC_totnod,LSC_totnod);


Bigfvec = zeros(numel,LSC_nnpe);

Cmat_els = zeros(numel,LSC_nnpe,LSC_nnpe);
% Dmat_els = zeros(numel,LSC_nnpe,LSC_nnpe);

% Total area of the domain
Area = 0;

% Matrices for interpolation
for iele = 1:1:numel
    
    % Shape function derivatives and determinant
    [dndx,dndy,detj] = LSC_sfder(iele, LSC_nnpe, LSC_nqpt, LSC_dN_qpt, LSC_conn, LSC_crds);
    
%     % Shape function derivatives and determinant
%     [dndx2,dndy2,detj2] = LSC_sfder2(iele, LSC_nnpe, LSC_nqpt, LSC_dN2_qpt, LSC_conn, LSC_crds);
    
    
    % Least Square mapping and pseudo force vector
    [Mmat, mvec, area] = LSC_mmat(LSC_nnpe, LSC_nqpt, LSC_N_qpt, LSC_wght, detj);
    
    % Overall area of the domain
    Area = Area + area;
    

    
    % Continuity constraint
    [Cmat,Cmat_norm] = LSC_cmat(LSC_nnpe, LSC_nqpt, LSC_wght, dndx,dndy,detj);
    
%     % Constant curvature gradient constraint
%     [Dmat, Dmat_norm, dvec] = LSC_dmat(LSC_nnpe, LSC_nqpt, LSC_wght, dndx,dndy,detj, div_kappa);
    
    
    % Store gradient operator
    Cmat_els(iele,1:LSC_nnpe,1:LSC_nnpe) = Cmat_norm;
    
    
%     % Store gradient operator
%     Dmat_els(iele,1:LSC_nnpe,1:LSC_nnpe) = Dmat_norm;
    
%     Nmat = Mmat + lambda*Cmat + beta*Dmat;
    Nmat = Mmat + lambda*Cmat ;
    
    
    % Store the pre-multiplier of the pseudo-force vector
    Bigfvec(iele,1:LSC_nnpe) = mvec; %+ dvec;
    
    % Assemble the local elemental constraints to the corresponding global
    % degrees of freedom
    [BigNmat] = LSC_assemblecons(iele,LSC_nnpe,LSC_conn,Nmat,BigNmat);
    
    [BigCmat] = LSC_assemblecons(iele,LSC_nnpe,LSC_conn,Cmat,BigCmat);
    
%     [BigDmat] = LSC_assemblecons(iele,LSC_nnpe,LSC_conn,Dmat,BigDmat);
    
    
end

% Interpolation functions at the center of the element
N0 = LSC_n0(LSC_eltyp);

clear k e


