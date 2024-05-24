function [Cmat, Cmat_norm, fvec] = cmat(nnpe, nqpt, wght, N_qpt, dndx, dndy, detj, w)
% LSC_cmat computes continuity constraint

Cmat = zeros(nnpe,nnpe);

fvec = zeros(nnpe,2);

% Area of the element
area = 0;

% Loop through the integration points
for i = 1:1:nqpt

    
    
    
    
    Nder = [ dndx(i,1:nnpe);
             dndy(i,1:nnpe); ];

     Cmat = Cmat + Nder'*Nder*detj(i)*wght(i);

     
     fvec = fvec + Nder'*detj(i)*wght(i);
     
     area = area + detj(i);

end

Cmat_norm = Cmat/area;

