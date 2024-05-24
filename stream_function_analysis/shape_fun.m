function [N] = shape_fun(k,e,eltyp)
% LSC_n computes shape functions

if eltyp == 3

    N = zeros(1,3);
    
    N(1) = 1-k-e;
    
    N(2) = k;
    
    N(3) = e;
    
elseif eltyp == 4
    
    N = zeros(1,4);
    
    N(1) = (1-k)*(1-e)/4;
   
    N(2) = (1+k)*(1-e)/4;
   
    N(3) = (1+k)*(1+e)/4;
   
    N(4) = (1-k)*(1+e)/4;

    
end



