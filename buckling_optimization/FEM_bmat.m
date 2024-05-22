function [B] = FEM_bmat(dN)

B = [ dN(1,1) 0 dN(1,2) 0  dN(1,3) 0;
    0 dN(2,1) 0 dN(2,2) 0 dN(2,3);
    dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3);];


