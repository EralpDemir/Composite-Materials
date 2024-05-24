% This script defines least squares and continuity constraint properties



%%%%%%%%%%%%%% Type of shape function that would be used for smoothing 
% 3: Linear Triangular shape functions
% 4: Linear Quadrilateral shape functions
% 8: Quadratic Quadrilateral shape functions
eltyp = 4;






% Ply to be analyzed
izPly = 1;


% Penalty factor for thickness gradient
alpha=0;

% Penalty factor for thickness divergence
beta=0;






% Tow width [m]
%w_tow=1/4*25.4e-3 *2; 
w_tow=1/4*25.4 * 2; 




