% This script defines the type of optimization problem based on in-plane or
% out-of-plane displacements






% Type of optimization
% Orthotropic materials only!!!
 

% % Options for optimization - CS
% options = setoptimoptions('algorithm','fminsearch','Display', 'iter',...
%     'GradObj','on', 'GradConstr', 'on', 'AlwaysHonorConstraints', 'all');
% 

% options = optimoptions(@fmincon,'Algorithm','sqp',...
%     'MaxIterations',1000, 'SpecifyObjectiveGradient',true, ...
%     'SpecifyConstraintGradient', true, 'MaxFunctionEvaluations', ...
%     40000,'ConstraintTolerance', 1e-8, 'StepTolerance', 1e-8,...
%     'OptimalityTolerance', 1e-8, 'HonorBounds',true, 'Display', 'iter');




% % Options for optimization - CS
% options = optimoptions(@fmincon,'Algorithm','sqp',...
%     'MaxIterations',1000, 'SpecifyObjectiveGradient',true, ...
%     'SpecifyConstraintGradient', true, 'MaxFunctionEvaluations', ...
%     40000,'ConstraintTolerance', 1e-8, 'StepTolerance', 1e-12,...
%     'OptimalityTolerance', 1e-8, 'HonorBounds',true, 'Display', 'iter');


% % Options for optimization - VS
% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%     'MaxIterations',500, 'SpecifyObjectiveGradient',true, ...
%     'SpecifyConstraintGradient', true, 'MaxFunctionEvaluations', ...
%     5000,'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-10,...
%     'OptimalityTolerance', 1e-10, 'Display', 'iter');

options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'MaxIterations',1000, 'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient', true, 'MaxFunctionEvaluations', ...
    5000,'ConstraintTolerance', 1e-8, 'StepTolerance', 1e-8,...
    'OptimalityTolerance', 1e-8, 'HonorBounds',true, 'Display', 'iter');

% Type of the stiffness of the design
% CON: constant stiffness
% VAR: variable stiffness
typSTIFF = 'CON';





