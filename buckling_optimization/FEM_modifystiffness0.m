function [K] = FEM_modifystiffness0(K, dofBCs, numBCs)
%FEM_applybcs modifies the stiffness matrix and force vector for the given BCs
% K = global coefficient matrix
% F = global force vector
% dofBCs = list of degrees of freedom with specified values
% valVals = specified values

for i = 1:1:numBCs
    
    id = dofBCs(i);
    
   
    K(id,:)=0;
    
    K(:,id)=0;
    
    
    K(id,id) =0;
    
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