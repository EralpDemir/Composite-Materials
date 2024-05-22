function [N] = LSC_n(k,e,eltyp)
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
   
elseif eltyp == 8
    
    N = zeros(1,8);
    
    N(5) = (1-k^2)*(1-e)/2;
    
    N(6) = (1+k)*(1-e^2)/2;
    
    N(7) = (1-k^2)*(1+e)/2;
    
    N(8) = (1-k)*(1-e^2)/2;
    
    N(1) = (1-k)*(1-e)/4 - (N(5)+N(8))/2;
    
    N(2) = (1+k)*(1-e)/4 - (N(5)+N(6))/2;
    
    N(3) = (1+k)*(1+e)/4 - (N(6)+N(7))/2;
    
    N(4) = (1-k)*(1+e)/4 - (N(7)+N(8))/2;
    
end



