
function [fe_v] = FEM_fipedlv(qL_v,nnpe,Fipedlv_ElNod,crds,conn,eltyp,dof,iele)

%FEM_fipedlv computes the element force for In-Plane Distributed Loads on the
% Edge of BFS rectangular plate elements along x-direction

% If force is zeros, assign all the force vectors zero!
fe_v = zeros(dof*nnpe,1);

switch eltyp
    
  
    case 'REC'
        
        % Dimensions of the element
        a = abs(crds(conn(iele,2),1)-crds(conn(iele,1),1))/2;
        b = abs(crds(conn(iele,4),2)-crds(conn(iele,1),2))/2;
        
        for iedge = 1:1:nnpe
            
            if Fipedlv_ElNod(iele,iedge) == 1
                
                if iedge == 1
                    % Edge load between the nodes 1-2
                    fe_v = fe_v + qL_v*a*[  0, 1, ...
                        0, 1, ...
                        0, 0, ...
                        0, 0, ]';
                    
                elseif iedge == 2
                    % Edge load between the nodes 2-3
                    fe_v = fe_v + qL_v*b*[  0, 0, ...
                        0, 1, ...
                        0, 1, ...
                        0, 0, ]';
                    
                elseif iedge == 3
                    % Edge load between the nodes 3-4
                    fe_v = fe_v +   qL_v*a*[0, 0, ...
                        0, 0, ...
                        0, 1, ...
                        0, 1, ]';
                    
                elseif iedge == 4
                    % Edge load between the nodes 4-1
                    fe_v = fe_v +   qL_v*b*[0, 1,  ...
                        0, 0, ...
                        0, 0,  ...
                        0, 1, ]';
                    
                end
                
            end
            
        end
        
    
        
        
end


