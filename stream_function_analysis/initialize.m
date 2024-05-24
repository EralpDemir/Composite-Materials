

%_____________ Find connectivity and node coordinates ___________________
switch eltyp
    
    case 3
        

        
        nnpe = 3;
        
        nqpt = 1;
        
        qpt = [ 1/3 1/3];
        
        wght = 1/2;
        
    case 4

        
        nnpe = 4;
        
        qpt = [ -1/sqrt(3) -1/sqrt(3);
            1/sqrt(3) -1/sqrt(3);
            1/sqrt(3)  1/sqrt(3);
            -1/sqrt(3)  1/sqrt(3)];
        
        nqpt = 4;
        
        wght = [ 1;
            1;
            1;
            1];
        
        
end

% Total number of nodes in the mesh
totnod = size(crds,1);

% Calculation of shape functions and their derivates at the quadrature points
N_qpt = zeros(nqpt,nnpe);
dN_qpt = zeros(nqpt,2,nnpe);

for iqpt = 1:1:nqpt
    
    k = qpt(iqpt,1);
    
    e = qpt(iqpt,2);
    
    N_qpt(iqpt,1:nnpe) = shape_fun(k,e,eltyp);
    
    dN_qpt(iqpt,1:2,1:nnpe) = shape_fun_der(k,e,eltyp);
    
end

N0=shape_fun(0,0,eltyp);



