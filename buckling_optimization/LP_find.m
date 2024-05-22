function [angPly, xi_A, xi_D] = LP_find(xi_As, xi_Ds, angPlys, ...
    noAngPoss, xi)
% LP_find finds the nearest possiblity in the LP space to the optimum value

distmin = 1e10;

for j = 1:1:noAngPoss
    
    % Lamination parameter values
    xi_A_(1:4,1) = xi_As(j,1:4)';
    
    xi_D_(1:4,1) = xi_Ds(j,1:4)';
    
    xis = [xi_A_(1); xi_A_(3); xi_D_(1); xi_D_(3) ];
    
    
%     if abs(xi_A_(2))<1e-10
%         if abs(xi_A_(4))<1e-10
%             if abs(xi_D_(2))<1e-10
%                if abs(xi_D_(4))<1e-10 
%             
    
%         dist = norm(xis(1:2) - xi(1:2));
%     dist = norm(xis(3:4) - xi(3:4));
    dist = norm(xis - xi);
    
    if dist < distmin
        distmin = dist;
        
        jmin = j;
        xi_A = xi_A_;
        xi_D = xi_D_;
        
    end
    
    
    
    
    
%                end
%             end
%         end
%     end
    
end

angPly(1,:) = angPlys(jmin,:);

% Error message
if distmin == 1e10
    error('COULD NOT FIND ANY NEIGHBOR! ');
end

return