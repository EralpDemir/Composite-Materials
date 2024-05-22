function [K] = FEM_modifystiffness(K, dofBCs, numBCs)
%FEM_applybcs modifies the stiffness matrix and force vector for the given BCs
% K = global coefficient matrix
% F = global force vector
% dofBCs = list of degrees of freedom with specified values
% valVals = specified values

for i = 1:1:numBCs
    
    id = dofBCs(i);
    

    
   K(id,id) = K(id,id) *1.0e+10;
    
end


%for i=1:1:numind
   
 %   if mod(i,dof) == 1
  %    K(i,:) =0;
   %   K(:,i) =0;
    %  K(i,i) =1;
  % end
       
   %if mod(i,dof) == 2
    %  K(i,:) =0;
    %  K(:,i) =0;
    %  K(i,i) =1; 
   %end
    
end