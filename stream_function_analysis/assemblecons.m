function [K] = assemblecons(iele,nnpe,conn,ke,K)

%LSC_assemblecons assembles the elemental stiffness matrices
%
% INPUTS:
% iele: number of the element
% nnpe: nodes per element
% conn: element connectivity
% ke: elemental stiffness matrix
% K: initial stiffness matrix (zero)
%
% OUTPUTS:
% K: assembled stiffness matrix on each element

for  i = 1:1:nnpe
    
    ii = conn(iele,i);
    
    for j = 1:1:nnpe
        
        jj = conn(iele,j);
        
        K(ii,jj)    =   K(ii,jj) +  ke(i,j);
        
    end
    
end
