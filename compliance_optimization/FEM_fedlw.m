
function [feL_w] = FEM_fedlw(qL,nnpe,Fedlw_ElNod,conn,crds,eltyp,dof,iele,numind)

%FEM_fedlw computes the element force for a distributed load on the
% edge of BFS rectangular plate elements

switch eltyp
    
    case 'TRI'
        
        feL_w = zeros(dof*nnpe,1);
        
        % Only one of the edges can have a distributed load
        for ied = 1:1:nnpe
            if Fedlw_ElNod(iele,ied) == 1
                error('Out-of-plane loads are NOT defined for TRI type of elements! Please check the definition for "Fsdlw_El"!')
            end
        end
        
    case 'HCT'
        
        feL_w = zeros(numind,1);
        
    case 'REC'
        
        feL_w = zeros(dof*nnpe,1);
        
        % Only one of the edges can have a distributed load
        for ied = 1:1:nnpe
            if Fedlw_ElNod(iele,ied) == 1
                error('Out-of-plane loads are NOT defined for REC type of elements! Please check the definition for "Fsdlw_El"!')
            end
        end
        
    case 'BFS'
        
        feL_w = zeros(dof*nnpe,1);
        
        % Dimensions of the element
        a = abs(crds(conn(iele,2),1) - crds(conn(iele,1),1))/2;
        b = abs(crds(conn(iele,4),2) - crds(conn(iele,1),2))/2;
        
        % Only one of the edges can have a distributed load
        for ied = 1:1:nnpe
            
            if Fedlw_ElNod(iele,ied) == 1
                
                if ied == 1
                    % Edge load between the nodes 1-2
                    feL_w = feL_w +  qL*[   a, a^2/3, 0, 0, ...
                        a, -a^2/3, 0, 0, ...
                        0, 0, 0, 0, ...
                        0, 0, 0, 0]';
                    
                elseif ied == 2
                    % Edge load between the nodes 2-3
                    feL_w =  feL_w +  qL*[  0, 0, 0, 0, ...
                        b, 0, b^2/3, 0, ...
                        b, 0, -b^2/3, 0, ...
                        0, 0, 0, 0]';
                    
                elseif ied == 3
                    % Edge load between the nodes 3-4
                    feL_w =  feL_w + qL*[ 0, 0, 0, 0, ...
                        0, 0, 0, 0, ...
                        a, -a^2/3, 0, 0, ...
                        a, a^2/3, 0, 0]';
                    
                elseif ied == 4
                    % Edge load between the nodes 4-1
                    feL_w =  feL_w + qL*[   b, 0, b^2/3, 0, ...
                        0, 0, 0, 0, ...
                        0, 0, 0, 0, ...
                        b, 0, -b^2/3, 0]';
                    
                end
                
            end
            
        end
        
end



