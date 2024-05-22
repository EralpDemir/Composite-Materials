% This script defines the type of optimization problem based on in-plane or
% out-of-plane displacements

% Type of optimization
% 1: optimize A only (in-plane displacement problems)
% 2: optimize D only (out-of-plane displacement problems)
typOPT = 2;


% Options for optimization
options = optimoptions('fmincon','Algorithm','interior-point',...
    'MaxIterations',10000, 'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, 'MaxFunctionEvaluations', ...
    10000,'ConstraintTolerance', 1e-9, 'StepTolerance', 2e-8,...
    'OptimalityTolerance', 1e-9, 'HonorBounds',true, 'Display', 'iter');


% Type of the stiffness of the design
% CON: constant stiffness
% VAR: variable stiffness
typSTIFF = 'VAR';


% Type of the laminate
% ORT: orthotropic laminate
% GEN: general laminate
typLAM = 'GEN';

