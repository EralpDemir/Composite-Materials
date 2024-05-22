
function [fe_u] = FEM_fipedlu(qL_u,nnpe,Fipedlu_ElNod,crds,conn,eltyp,dof,iele,numind)

%FEM_fipedlu computes the element force for In-Plane Distributed Loads on the
% Edge of BFS rectangular plate elements along x-direction

% If force is zeros, assign all the force vectors zero!
fe_u = zeros(dof*nnpe,1);

switch eltyp
    
    case 'TRI'
        
        for iedge = 1:1:nnpe
            
            if Fipedlu_ElNod(iele,iedge) == 1
                
                x1 = crds(conn(iele,1),1);  y1 = crds(conn(iele,1),2);
                
                x2 = crds(conn(iele,2),1);  y2 = crds(conn(iele,2),2);
                
                x3 = crds(conn(iele,3),1);  y3 = crds(conn(iele,3),2);
                
                if iedge == 1
                    
                    % Edge load between the nodes 1-2
                    
                    L = sqrt((x2-x1)^2+(y2-y1)^2);
                    
                    fe_u = fe_u + qL_u*L*[  0.5, 0, ...
                        0.5, 0, ...
                        0, 0,]';
                    
                elseif iedge == 2
                    % Edge load between the nodes 2-3
                    
                    L=sqrt((x2-x3)^2+(y2-y3)^2);
                    
                    fe_u = fe_u + qL_u*L*[  0, 0, ...
                        0.5, 0, ...
                        0.5, 0,]';
                    
                elseif iedge == 3
                    % Edge load between the nodes 3-1
                    
                    L = sqrt((x3-x1)^2+(y3-y1)^2);
                    
                    fe_u = fe_u + qL_u*L*[  0.5, 0, ...
                        0, 0, ...
                        0.5, 0 ]';
                end
                
            end
            
        end
        
    case 'HCT'
        
        fe_u = zeros(numind,1);
        
        for inod = 1:1:nnpe
            
            if Fipedlu_ElNod(iele,inod) == 1
                error('In-plane loads are NOT defined for HCT type of elements! Please check the definition for "Fipedlu_ElNod"!')
            end
            
        end
        
    case 'REC'
        
        % Dimensions of the element
        a = abs(crds(conn(iele,2),1)-crds(conn(iele,1),1))/2;
        b = abs(crds(conn(iele,4),2)-crds(conn(iele,1),2))/2;
        
        for iedge = 1:1:nnpe
            
            if Fipedlu_ElNod(iele,iedge) == 1
                
                if iedge == 1
                    % Edge load between the nodes 1-2
                    fe_u = fe_u + qL_u*a*[  1, 0, ...
                        1, 0, ...
                        0, 0, ...
                        0, 0, ]';
                    
                elseif iedge == 2
                    % Edge load between the nodes 2-3
                    fe_u = fe_u + qL_u*b*[  0, 0, ...
                        1, 0, ...
                        1, 0, ...
                        0, 0, ]';
                    
                elseif iedge == 3
                    % Edge load between the nodes 3-4
                    fe_u = fe_u + qL_u*a*[0, 0, ...
                        0, 0, ...
                        1, 0, ...
                        1, 0, ]';
                    
                elseif iedge == 4
                    % Edge load between the nodes 4-1
                    fe_u = fe_u +   qL_u*b*[1, 0,  ...
                        0, 0, ...
                        0, 0,  ...
                        1, 0, ]';
                    
                end
                
            end
            
        end
        
    case 'BFS'
        
        for iedge = 1:1:nnpe
            
            if Fipedlu_ElNod(iele,iedge) == 1
                error('In-plane loads are NOT defined in BFS! Please check the definition for "Fipedlu_ElNod"')
            end
            
        end
        
end


