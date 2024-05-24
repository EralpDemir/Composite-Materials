function [Nmat, Cmat, Dmat, fvec, area] = nmat(nnpe, nqpt, wght, N_qpt, dndx, dndy, detj, angs)
% LSC_cmat computes continuity constraint

Nmat = zeros(nnpe,nnpe);

Cmat = zeros(nnpe,nnpe);


Dmat = zeros(nnpe,nnpe);

fvec = zeros(nnpe,1);

% Area of the element
area = 0;


% % angle at the center of the element
% ang0 = N0 * angs;

% Loop through the integration points
for i = 1:1:nqpt

    N(1,1:nnpe)=N_qpt(i,:);
    
    % Angle at the quadrature point
    ang = N * angs;
    
    
    Nder = [ dndx(i,1:nnpe);
             dndy(i,1:nnpe); ];
    

    % Gradient of angles     
    dang = Nder * angs;
    
    
    

    
    n_bar = [-sin(ang), cos(ang)];
    s_bar = [cos(ang), sin(ang)];


    
    
    
    % interpolation matrix
%      Nmat = Nmat + N'*cN*Nder*detj(i)*wght(i);
     
    Nmat = Nmat + (s_bar*Nder)' * (s_bar*Nder) * detj(i) * wght(i);
    
    Cmat = Cmat  + Nder'  * Nder * detj(i) * wght(i);
    
    Dmat = Dmat  + Nder'  * [1, 1]' * [1, 1] * Nder * detj(i) * wght(i);

%      Nmat = Nmat + Nder'*cN'*cN*Nder*detj(i)*wght(i);
     
%     Nmat = Nmat + Nder'*[1, 1]'*cN*Nder*detj(i)*wght(i);
    
    % pseudo-force vector
%     fvec = fvec + N'*detj(i)*wght(i)*cf*dang;
    
    fvec = fvec + (s_bar*Nder)' * detj(i) * wght(i) * n_bar * dang;
    
    area = area + detj(i);

end



