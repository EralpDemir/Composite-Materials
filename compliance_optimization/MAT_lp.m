function [xi_A,xi_B,xi_D] = MAT_lp(nPly, angPly, zPly, t)

%MAT_LP     calculates the lamination parameters space 
% (vector of harmonics)
% 
% INPUTS:
% nPly: total number of laminate plies
% angPly: fiber angle at each lamina layer
% zPly: start and end of each lamina layer 
% t: total thickness of laminate
% 
% OUTPUTS:
% xi_A: parameter space for Lamination Parameter-A
% xi_B: parameter space for Lamination Parameter-B
% xi_D: parameter space for Lamination Parameter-D
%

conv = pi/180; % Conversion to radians

% Lamination parameters for different thicknesses
xi_A = zeros(4,1); 
xi_B = zeros(4,1); 
xi_D = zeros(4,1);
for i = 1:1:nPly
    
    theta = angPly(i)*conv; 
    
    z = zPly(i,2); % upper surface
    z_= zPly(i,1); % lower surface
    
    % Calculation of Lamination Parameters (SUM)
    xi_A = xi_A + ...
        [cos(2*theta); sin(2*theta); cos(4*theta); sin(4*theta)]...
        *(z-z_)/t;

    xi_B = xi_B + ...
        [cos(2*theta); sin(2*theta); cos(4*theta); sin(4*theta)]...
        *2*(z^2-z_^2)/t^2;
    
    xi_D = xi_D + ...
        [cos(2*theta); sin(2*theta); cos(4*theta); sin(4*theta)]...
        *4*(z^3-z_^3)/t^3;
end

