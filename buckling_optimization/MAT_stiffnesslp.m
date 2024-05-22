% Eralp Demir
% Sept. 29th, 2016

function [A_vec, D_vec] = MAT_stiffnesslp(xi_A, xi_D, t, U)



% A matrix components - lamina
A_vec = t*U*[1; xi_A];



% D matrix components - lamina
D_vec = t^3/12*U*[1; xi_D];



