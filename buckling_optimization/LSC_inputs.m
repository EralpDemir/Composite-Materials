% This script defines least squares and continuity constraint properties

%%%%%%%%%%%%%% Constraint penalty parameters (continutity constraint)
% Curvature magnitude constraint
alpha = 100;
% alpha = 2200; % 1/400
% alpha = 3700; % 1/600

% % Curvature gradient constraint
% beta = 0.01;
% 
% % Curvature gradient value
% div_kappa=0;

%%%%%%%%%%%%%% Type of shape function that would be used for smoothing 
% 3: Linear Triangular shape functions
% 4: Linear Quadrilateral shape functionse
% 8: Quadratic Quadrilateral shape functions
LSC_eltyp = 4;


