function [A_vec, D_vec] = OPT_tangent(t, U)

%OPT_Tangent    calculates analytical 1st order sensitivities
% 
% INPUTS:
% t: total thickness of laminate
% U: stiffness lamination parameters matrix
% U_: shear lamination parameters matrix
%
% OUTPUTS:
% A_vec: 1st order sensitivity A matrix
% B_vec: 1st order sensitivity B matrix
% D_vec: 1st order sensitivity D matrix
% H_vec: 1st order sensitivity H matrix

% A matrix components - lamina
A_vec(:,1) = t*U*[0 1 0 0 0]';
A_vec(:,2) = t*U*[0 0 1 0 0]';
A_vec(:,3) = t*U*[0 0 0 1 0]';
A_vec(:,4) = t*U*[0 0 0 0 1]';

% % B matrix components - lamina
% B_vec(:,1) = t^2/4*U*[0 1 0 0 0]';
% B_vec(:,2) = t^2/4*U*[0 0 1 0 0]';
% B_vec(:,3) = t^2/4*U*[0 0 0 1 0]';
% B_vec(:,4) = t^2/4*U*[0 0 0 0 1]';

% D matrix components - lamina
D_vec(:,1) = t^3/12*U*[0 1 0 0 0]';
D_vec(:,2) = t^3/12*U*[0 0 1 0 0]';
D_vec(:,3) = t^3/12*U*[0 0 0 1 0]';
D_vec(:,4) = t^3/12*U*[0 0 0 0 1]';

% % H matrix components - lamina
% H_vec(:,1) = t*U_*[0 1 0 0 0]';
% H_vec(:,2) = t*U_*[0 0 1 0 0]';
% H_vec(:,3) = t*U_*[0 0 0 1 0]';
% H_vec(:,4) = t*U_*[0 0 0 0 1]';

