% This script plots the lamination parameters and angle distribution for REC and BFS

if typOPT==1
    
    for i = 1:4
        OPT_outputlp(nx, ny, L_x, ratio, All_xi_A(i,:))
        caxis([-1 1]);
        axis equal
    end
    
elseif typOPT == 2
    
    for i = 1:4
        OPT_outputlp(nx, ny, L_x, ratio, All_xi_D0(i,:))
        caxis([-1 1]);
        axis equal
    end
    
end

% Plot angle distribution
for i = 1:nPly/2
    OPT_outputel(conn,crds,numel,nnpe,angPly,i)
    fnam_var = ['OptAng-L',num2str(i)];
    fname = [Case_Name,fnam_var,fnam_Format];
end
