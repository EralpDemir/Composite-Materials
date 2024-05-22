% Eralp Demir
% Sept. 29th, 2016

function [A_vec, B_vec, D_vec] = MAT_stiffness(nPly, angPly, zPly, t, U)

% The vector of harmonics
[xi_A, xi_B, xi_D] = MAT_lp(nPly, angPly, zPly, t);

% A matrix components - lamina
A_vec = t*U*[1; xi_A];

% B matrix components - lamina
B_vec = t^2/4*U*[0; xi_B];

% D matrix components - lamina
D_vec = t^3/12*U*[1; xi_D];



