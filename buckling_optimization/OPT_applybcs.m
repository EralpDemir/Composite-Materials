function [F] = OPT_applybcs(F, dofBCs, valBCs, numBCs)
%FEM_applybcs modifies the stiffness matrix and force vector for the given BCs
% K = global coefficient matrix
% F = global force vector
% dofBCs = list of degrees of freedom with specified values
% valVals = specified values

for i = 1:1:numBCs
    
    id = dofBCs(i);
    

    F(id)    = valBCs(i);
    
end