function [Mmat, mvec, area] = LSC_mmat(nnpe, nqpt, N_qpt, LSC_wght, detj)
% LSC_mmat computes least Square mapping and pseudo force vector

Mmat = zeros(nnpe,nnpe);

mvec = zeros(nnpe,1);

% total area of the element
area = 0;

% Loop through the integration points
for i=1:1:nqpt
       
    N(1,1:nnpe)=N_qpt(i,:);

    Mmat = Mmat + N'*N*detj(i)*LSC_wght(i);

    mvec = mvec + N'*detj(i)*LSC_wght(i);
    
    area = area + detj(i);

end

