function [N0] = LSC_n0(LSC_eltyp)
% LSC_n0 finds interpolation functions at the center of the element

% Center coordinate of the element
if LSC_eltyp == 3
    k = 1/3;
    e = 1/3;
elseif LSC_eltyp == 4
    k = 0;
    e = 0;
elseif LSC_eltyp == 8
    k = 0;
    e = 0;
end

% shape functions
N0 = LSC_n(k,e,LSC_eltyp);

