
function [feS] = FEM_fsdlw(q,Fsdlw_El,crds,conn,eltyp,nnpe,dof,iele,numind,nqpts)

%FEM_fsdlw This function computes the element force for a distributed load

switch eltyp
    
    case 'TRI'
        
        feS = zeros(dof*nnpe,1);
        if Fsdlw_El(iele) == 1
            error('Out-of-plane loads are NOT defined for TRI type of elements! Please check the definition for "Fsdlw_El"!')
        end
        
    case 'HCT'
        
        if Fsdlw_El(iele) == 1
            
            % Global node numbers for each element
            elnodes = conn(iele,:);
            
            % Node coordinates for each element global nodes regarding connectivity
            elcrds(:,:) = crds(elnodes,:);
            
            % First 3 coordinates are the vertex coordinates for each element
            % Note that these are local nodes in each element. global node 3 1 2 6 4 5
            % is related to local node 1 2 3 4 5 6
            
            x1 = elcrds(1,1);
            y1 = elcrds(1,2);
            x2 = elcrds(2,1);
            y2 = elcrds(2,2);
            x3 = elcrds(3,1);
            y3 = elcrds(3,2);
            
            % Compute quad point locations
            [xi, eta, w] = FEM_quad_trgl(nqpts);
            
            % compute the area of the triangle
            d12x = x1-x2; d12y = y1-y2;
            d31x = x3-x1; d31y = y3-y1;
            area = 0.5*(d31x*d12y - d31y*d12x); % area of each element
            
            % interior node
            x7 = (x1+x2+x3)/3.0;
            y7 = (y1+y2+y3)/3.0;
            
            % initialize the element bending stiffness matrix and right-hand side
            feS = zeros(numind,1);
            
            % Deformation modes
            am = zeros(30,numind);
            for mode = 1:numind
                
                dofs = zeros(1,numind);
                
                dofs(mode) = 1.0;
                
                [a] = FEM_hctsys (x1,y1, x2,y2, x3,y3, dofs);
                
                for i = 1:30
                    am(i,mode) = a(i);
                end
                
            end
            
            % Sub-triangles
            for isub = 1:3
                
                % Quad points of each sub-element
                for iquad = 1:nqpts   % numerical quadrature
                    
                    psi = zeros(1,numind);
                    
                    % For each deformation mode
                    for imode = 1:numind
                        
                        c(:,1) = am(:,imode);
                        
                        if isub == 1         % triangle 1-2-7
                            
                            x = x1 + (x2-x1)*xi(iquad) + (x7-x1)*eta(iquad);
                            y = y1 + (y2-y1)*xi(iquad) + (y7-y1)*eta(iquad);
                            
                            psi(imode) = c(1) + c(2)*x   + c(3)*y ...
                                + c(4)*x^2 + c(5)*x*y   + c(6)*y^2 ...
                                + c(7)*x^3 + c(8)*x^2*y + c(9)*x*y^2 + c(10)*y^3;
                            
                        elseif isub == 2   % triangle 2-3-7
                            
                            x = x2 + (x3-x2)*xi(iquad) + (x7-x2)*eta(iquad);
                            y = y2 + (y3-y2)*xi(iquad) + (y7-y2)*eta(iquad);
                            
                            psi(imode) = c(11) + c(12)*x   +c(13)*y ...
                                + c(14)*x^2 + c(15)*x*y   + c(16)*y^2 ...
                                + c(17)*x^3 + c(18)*x^2*y + c(19)*x*y^2 + c(20)*y^3;
                            
                        else            % triangle 3-1-7
                            
                            x = x3 + (x1-x3)*xi(iquad) + (x7-x3)*eta(iquad);
                            y = y3 + (y1-y3)*xi(iquad) + (y7-y3)*eta(iquad);
                            
                            psi(imode) = c(21) + c(22)*x + c(23)*y ...
                                + c(24)*x^2 + c(25)*x*y   + c(26)*y^2 ...
                                + c(27)*x^3 + c(28)*x^2*y + c(29)*x*y^2 + c(30)*y^3;
                            
                        end
                        
                    end
                    
                    % Area and weight of integraion of each subtriangle
                    cf = area*w(iquad)/3.0;
                    
                    for k = 1:numind
                        
                        feS(k) = feS(k) + q*psi(k)*cf;
                        
                    end
                    
                end
                
            end
            
        end
        
    case 'REC'
        
        feS = zeros(dof*nnpe,1);
        
        if Fsdlw_El(iele) == 1
            error('Out-of-plane loads are NOT defined for REC type of elements! Please check the definition for "Fsdlw_El"!')
        end
        
    case 'BFS'
        
        feS = zeros(dof*nnpe,1);
        
        % Dimensions of the element
        a = abs(crds(conn(iele,2),1)-crds(conn(iele,1),1))/2;
        b = abs(crds(conn(iele,4),2)-crds(conn(iele,1),2))/2;
        
        if Fsdlw_El(iele) == 1
            
            % Input force vector
            feS = a*b*q/9*[9, 3*a, 3*b, a*b,  9, -3*a, 3*b, -(a*b), ...
                9, -3*a, -3*b, a*b, 9, 3*a, -3*b, -(a*b)]';
            
        end
end

end
