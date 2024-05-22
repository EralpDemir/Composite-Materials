function  [x_c] = LSC_LP0(numel, nnpe, conn, x_nodes, N0)
%LSC_angel0 calculate the angles at the element centers

x_c = zeros(numel,4);

for i = 1:1:4
    
    ang_vals(:,1) = x_nodes(:,i);
    
    for iele = 1:1:numel

        ang_vals_ele = zeros(nnpe,1);
        
        for inod = 1:1:nnpe

            NodeNo = conn(iele,inod);

            ang_vals_ele(inod,1) = ang_vals(NodeNo);

        end

        x_c(iele,i) = N0*ang_vals_ele;

    end
        
end
    
