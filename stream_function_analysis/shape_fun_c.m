function [N0, dN0] = shape_fun_c(eltyp)
% LSC_n0 finds interpolation functions at the center of the element

% Center coordinate of the element
if eltyp == 3
    
    k = 1/3;
    e = 1/3;
    
elseif eltyp == 4
    
    k = 0;
    e = 0;
    
end

% shape functions
N0 = shape_fun(k,e,eltyp);

% shape functions
dN0 = shape_fun_der(k,e,eltyp);



