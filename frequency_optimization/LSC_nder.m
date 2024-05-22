function [Nder] = LSC_nder(k,e,LSC_eltyp)

if LSC_eltyp == 3

    Nder(1,1) = -1;
    
    Nder(1,2) = 1;
    
    Nder(1,3) = 0;
    
    Nder(2,1) = -1;
    
    Nder(2,2) = 0;
    
    Nder(2,3) = 1;

elseif LSC_eltyp==4
    
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
    
elseif LSC_eltyp == 8
    
    Nder(1,5) = (-k)*(1-e);
    %
    Nder(1,6) = (1)*(1-e^2)/2;
    %
    Nder(1,7) = (-k)*(1+e);
    %
    Nder(1,8) = (-1)*(1-e^2)/2;
    
    Nder(1,1) = (-1)*(1-e)/4 - (Nder(1,8) + Nder(1,5))/2;
    %
    Nder(1,2) = (1)*(1-e)/4 - (Nder(1,5) + Nder(1,6))/2;
    %
    Nder(1,3) = (1)*(1+e)/4 - (Nder(1,6) + Nder(1,7))/2;
    %
    Nder(1,4) = (-1)*(1+e)/4 - (Nder(1,7) + Nder(1,8))/2;
    %
    Nder(2,5) = (1-k^2)*(-1)/2;
    %
    Nder(2,6) = (1+k)*(-e);
    %
    Nder(2,7) = (1+k^2)*(1)/2;
    %
    Nder(2,8) = (1-k)*(-e);
    
    Nder(2,1) = (1-k)*(-1)/4 - (Nder(2,8) + Nder(2,5))/2;
    %
    Nder(2,2) = (1+k)*(-1)/4 - (Nder(2,5) + Nder(2,6))/2;
    %
    Nder(2,3) = (1+k)*(1)/4 - (Nder(2,6) + Nder(2,7))/2;
    %
    Nder(2,4) = (1-k)*(1)/4 - (Nder(2,7) + Nder(2,8))/2;
    
end


