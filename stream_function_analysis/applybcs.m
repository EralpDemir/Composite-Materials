function [K, F] = applybcs(K, F, dofBCs, valBCs, numBCs)
%FEM_applybcs modifies the stiffness matrix and force vector for the given BCs
% K = global coefficient matrix
% F = global force vector
% dofBCs = list of degrees of freedom with specified values
% valVals = specified values

for i = 1:1:numBCs
    
    id = dofBCs(i);
    
    K(id,id) = K(id,id) *1.0e+10;
    F(id)    = valBCs(i)*K(id,id);
    
end