
function [f] = FEM_fpl(f, dof, Fpl_NodComp, p)

% FEM_fpl computes the element force for a distributed load on BFS
% rectangular plate elements
%
% f: global force vector
% b: dimension of rectangular element along y-direction
% q: distributed load over the area

if not(p==0)
    for i = 1:1:size(Fpl_NodComp,1)
        
        numnod = Fpl_NodComp(i,1);
        comp = Fpl_NodComp(i,2);
        
        f(dof*(numnod-1)+comp) = p;
        
    end
    
end

end
