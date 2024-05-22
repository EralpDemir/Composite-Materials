function [xi_As, xi_Bs, xi_Ds, noAngPoss, angPlys] = ...
    LP_space(ang_range, t, Ltype, nPly, zPly)

%LP_Space   calculates the lamination parameters space
%
% INPUTS:
% ang_range: range of possible fiber angles in each lamina
% nPlys: total number of plies
% t: total thickness of laminate
% Ltype: type of laminate
%
% OUTPUTS:
% xi_As: parameter space for Lamination Parameter-A
% xi_Bs: parameter space for Lamination Parameter-B
% xi_Ds: parameter space for Lamination Parameter-D
% noPlyPoss: possible number of plies
% noAngPoss: possible number of angles
% zPlys: z coordinate of plies
% angPlys: orientation of fibers in plies
%

%---------------------------------------------------------------------
%%%%%%%% uncomment to generate the figure lamination space
% www = 14.3; % for ratio 1
% hhh = 12.6; % for ratio 1
% figure('Units','centimeters','Position',[0 0 www hhh])
% axis off
% axis equal
% set(gcf,'paperunits','centimeters','papersize',[www hhh],...
%     'paperposition',[0 0 www hhh])
% positionVector1 = [2, 1.8, 12, 10.5]; % for ratio 1
% axes('units','centimeters','position',positionVector1);
% xticks(-1:0.2:1); yticks(-1:0.2:1)
% xlabel(['\color{blue}\xi_2^A',' \color{black}, ','\color{red}\xi_2^B',...
%     ' \color{black}, ','\color{black}\xi_2^D'],'FontName','Cambria Math',...
%     'FontWeight','bold','FontSize',16,'interpreter','tex');
% ylabel(['\color{blue}\xi_4^A',' \color{black}, ','\color{red}\xi_4^B',...
%     ' \color{black}, ','\color{black}\xi_4^D'],'FontName','Cambria Math',...
%     'FontWeight','bold','FontSize',16,'interpreter','tex');
% grid on
% hold on
%---------------------------------------------------------------------

zPly = [zPly(1,:,1); zPly(1,:,2)]';
angPlys = [];

%___________________________ Single Layer Laminate _______________________
if Ltype == 1
    
    % Form all the combinations of angles

    ANG = LP_permn(ang_range, 2*nPly);
    if size(ANG,1) == 1
        ANG = ANG';
    end
    
    %___________________________ Symmetric Laminate ______________________
elseif Ltype == 2
    

    ANG = LP_permn(ang_range, nPly);
    if size(ANG,1) == 1
        ANG = ANG';
    end
    
    % Re-arrange the possibilities and add the SYMMETRIC Layup
    ANG = [ANG ANG];
    
    %_________________________ Balanced Laminate ________________________
elseif Ltype == 3
    
    % Form all the combinations of angles

    ANG = LP_permn(ang_range, 2*nPly);
    if size(ANG,1) == 1
        ANG = ANG';
    end
    ANG = LP_permute_blnc(ANG, nPly);
    
    %_______________________ Symmetric & Balanced Laminate ______________
elseif Ltype == 4
    
    if rem(nPly/2,2) == 0
        
        % Form all the combinations of angles

        ANG = LP_permn(ang_range, nPly/2);
        if size(ANG,1) == 1
            ANG = ANG';
        end
        
        ANG_neg = -ANG;
        ANG_zero = ANG_neg == 0;
        ANG_ninty = ANG_neg == 90;
        ANG_neg(ANG_zero) = -90;
        ANG_neg(ANG_ninty) = 0;
        ANG_comb = [ANG, ANG_neg];
        P = randperm(nPly/2,nPly/2);
        clear ANG;
        % Re-arrange the possibilities
        for i = 1:nPly/2
            ANG(:,i) = ANG_comb(:,P(i));
        end
        
        ANG = [ANG -ANG];
        
    else
        error('defined # of plies (nPly) is not proper for this case!');
    end
    
    %___________________________ Asymmetric Laminate _____________________
elseif Ltype == 5
    
    % Form all the combinations of angles

    ANG = LP_permn(ang_range, nPly);
    if size(ANG,1) == 1
        ANG = ANG';
    end
    
    % Re-arrange the possibilities and add the SYMMETRIC Layup
    ANG = [ANG -ANG];
    
    %___________________________ General Laminate ___________________________
elseif Ltype == 6
    
    % Form all the combinations of angles

    ANG = LP_permn(ang_range, nPly); 
    if size(ANG,1) == 1
        ANG = ANG';
    end
    
     %___________________________ Specially Orthotropic Laminate ___________________________
     % A = [q / -q / -q / q], symmetric & balanced
     % [A / -A], anti-symmetric
elseif Ltype == 7
    
    
    ANG = ang_range';
    
    ANG = [-ANG, ANG, ANG, -ANG, ANG, -ANG, -ANG, ANG];
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%% Calculation and Plotting LPs %%%%%%%%%%%%%%%%%%%
noAngPoss_ = size(ANG,1);

angPlys_(1:size(ANG,1),1:size(ANG,2)) = ANG;

j=0;

for i = 1:1:noAngPoss_
    
    
    
    angPly = ANG(i,:);
    
    % Calculate lamination parameters
    [xi_A, xi_B, xi_D] = MAT_lp(nPly, angPly, zPly, t);
    
    x=[xi_A(1); xi_A(3); xi_D(1); xi_D(3)];
    
    [c]= OPT_constraint_fun(x);
    
    if all(c<=0)
    
        j=j+1;
        
        
        %             % Store and plot lamination parameter alternatives
        xi_As(j,:) = xi_A;
        %           plot(xi_A(2),xi_A(4),'b.','MarkerSize', 10)

        xi_Bs(j,:) = xi_B;
        %           plot(xi_B(2),xi_B(4),'r.','MarkerSize', 22)

        xi_Ds(j,:) = xi_D;
        %          plot(xi_D(2),xi_D(4),'k.','MarkerSize', 10)

        xis_(j,1:12) = [xi_A;xi_B;xi_D];
        
        angPlys(j,1:size(angPly,2))=angPly;
    
    end


end
%     saveas(gcf,'LPspace13.pdf')
%     print('LPspace24','-dpdf','-r1200')

noAngPoss = j;

return
end
