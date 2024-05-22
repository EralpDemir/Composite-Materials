% Nov. 21, 2016
% Lamination Parameter Space

function [xi_A,xi_B,xi_D] = MAT_lp(nPly, angPly, zPly, t)

conv=pi/180;

% nPly: total number of lamina layers
% angPly: fiber angle at each lamina layer
% zPly: start and end of each lamina layer (nPly x 2)


% Lamination parameters for different thicknesses
xi_A=zeros(4,1); 
xi_B=zeros(4,1); 
xi_D=zeros(4,1);
%
for i=1:1:nPly
    
    % Conversion to radians
    theta=angPly(i)*conv;
    
    % Thickness of lamina layer
    z = zPly(i,2);
    z_= zPly(i,1);
    
    % Calculation of Lamination Parameters (SUM)

    % Parameter space for Lamination Parameter-A
    xi_A = xi_A + ...
        [cos(2*theta);sin(2*theta);cos(4*theta);sin(4*theta)]...
        /t*(z-z_);

    % Parameter space for Lamination Parameter-B
    xi_B = xi_B + ...
        [cos(2*theta);sin(2*theta);cos(4*theta);sin(4*theta)]...
        *2/t^2*(z^2-z_^2);
    
    % Parameter space for Lamination Parameter-C
    xi_D = xi_D + ...
        [cos(2*theta);sin(2*theta);cos(4*theta);sin(4*theta)]...
        *4/t^3*(z^3-z_^3);
    
end


