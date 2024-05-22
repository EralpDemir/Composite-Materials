function [Cmat, Cmat_norm] = LSC_cmat(LSC_nnpe, LSC_nqpt, LSC_wght, dndx, dndy, detj)
% LSC_cmat computes continuity constraint

Cmat = zeros(LSC_nnpe,LSC_nnpe);

% Area of the element
area = 0;

% Loop through the integration points
for i = 1:1:LSC_nqpt

    Nder = [ dndx(i,1:LSC_nnpe);
             dndy(i,1:LSC_nnpe); ];

     Cmat = Cmat + Nder'*Nder*detj(i)*LSC_wght(i);

     area = area + detj(i);

end

Cmat_norm = Cmat/area;

