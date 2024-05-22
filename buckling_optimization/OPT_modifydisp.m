function [d] = OPT_modifydisp(d, dofBCs, numBCs)
%FEM_applybcs modifies the stiffness matrix and force vector for the given BCs
% K = global coefficient matrix
% F = global force vector
% dofBCs = list of degrees of freedom with specified values
% valVals = specified values

for i = 1:1:numBCs
    
    id = dofBCs(i);
    
   
    d(id)=0;
    
end



   
    
return