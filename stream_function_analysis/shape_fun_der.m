function [Nder] = shape_fun_der(k,e,eltyp)

if eltyp == 3

    Nder(1,1) = -1;
    
    Nder(1,2) = 1;
    
    Nder(1,3) = 0;
    
    Nder(2,1) = -1;
    
    Nder(2,2) = 0;
    
    Nder(2,3) = 1;

elseif eltyp==4
    
    Nder(1,1) = (-1)*(1-e)/4;
    %
    Nder(1,2) = (1)*(1-e)/4;
    %
    Nder(1,3) = (1)*(1+e)/4;
    %
    Nder(1,4) = (-1)*(1+e)/4;
    %
    Nder(2,1) = (1-k)*(-1)/4;
    %
    Nder(2,2) = (1+k)*(-1)/4;
    %
    Nder(2,3) = (1+k)*(1)/4;
    %
    Nder(2,4) = (1-k)*(1)/4;
    

    
end


