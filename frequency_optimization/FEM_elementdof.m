function [index] = FEM_elementdof(nd,nnpe,dof,eltyp)
%----------------------------------------------------------
%Compute system dofs associated with each element
%
%  Variable Description:
%     index - system dof vector associated with element "iel"
%     iel - element number whose system dofs are to be determined
%     nnpe - number of nodes per element
%     dof - number of dofs per node
%-----------------------------------------------------------
switch eltyp
    
    case 'TRI'
        index = [];
        for i = 1:nnpe
            index = [index (nd(i)*dof)-1 nd(i)*dof];
        end
        
    case 'REC'
        index = [];
        for i = 1:nnpe
            index = [index (nd(i)*dof)-1 nd(i)*dof];
        end
        
    case 'BFS'
        k=0;
        for i=1:nnpe
            start = (nd(i)-1)*dof;
            for j=1:dof
                k=k+1;
                index(k)=start+j;
            end
        end
        
end