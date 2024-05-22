function [lpsi] = FEM_ElStrDispB(iele, crds,conn, numind)

% This function calculates HCT triangular plate element stiffness

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

% interior node
x7=(x1+x2+x3)/3.0; y7=(y1+y2+y3)/3.0;

% Deformation modes
for mode=1:numind
    
    dofs = zeros(1,numind);
    
    dofs(mode) = 1.0;
    
    [a] = FEM_hctsys (x1,y1, x2,y2, x3,y3, dofs);
    
    for i=1:30
        am(i,mode) = a(i);
    end
    
end

% Sub-triangles
for isub = 1:3
    
    psi = zeros(1,numind);
    
    lpsi = zeros(3,numind);
    
    % For each deformation mode
    for imode = 1:numind
        
        c(:,1) = am(:,imode);
        
        if isub == 1        % triangle 1-2-7
            
            x = x7;
            y = y7;
            
            psi(imode) = c(1) + c(2)*x  + c(3)*y ...
                + c(4)*x^2 + c(5)*x*y   + c(6)*y^2 ...
                + c(7)*x^3 + c(8)*x^2*y + c(9)*x*y^2 + c(10)*y^3;
            
            lpsi(1:3,imode) = [   c(4)*2 + c(7)*6*x + c(8)*2*y;
                c(6)*2 + c(9)*2*x + c(10)*6*y;
                c(5) + c(8)*2*x + c(9)*2*y];
            
        elseif isub == 2   % triangle 2-3-7
            
            x = x7;
            y = y7;
            
            psi(imode) = c(11) + c(12)*x   +c(13)*y ...
                + c(14)*x^2 + c(15)*x*y   + c(16)*y^2 ...
                + c(17)*x^3 + c(18)*x^2*y + c(19)*x*y^2 + c(20)*y^3;
            
            lpsi(1:3,imode) = [   c(14)*2 + c(17)*6*x + c(18)*2*y;
                c(16)*2 + c(19)*2*x + c(20)*6*y;
                c(15) + c(18)*2*x + c(19)*2*y];
            
        else            % triangle 3-1-7
            
            x = x7;
            y = y7;
            
            psi(imode) = c(21) + c(22)*x + c(23)*y ...
                + c(24)*x^2 + c(25)*x*y   + c(26)*y^2 ...
                + c(27)*x^3 + c(28)*x^2*y + c(29)*x*y^2 + c(30)*y^3;
            
            lpsi(1:3,imode) = [ c(24)*2 + c(27)*6*x + c(28)*2*y;
                c(26)*2 + c(29)*2*x + c(30)*6*y;
                c(25) + c(28)*2*x + c(29)*2*y];
            
        end
    end
    
end

return