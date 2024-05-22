function  [lp_c] = LSC_lpel0(numel, nnpe, conn, lpNodes, N0)
%LSC_angel0 calculate the angles at the element centers

lp_c = zeros(1,numel);
    
lp_vals(:,1) = lpNodes(:);

for iele = 1:1:numel

    lp_vals_ele = zeros(nnpe,1);

    for inod = 1:1:nnpe

        NodeNo = conn(iele,inod);

        lp_vals_ele(inod,1) = lp_vals(NodeNo);

    end

    lp_c(1,iele) = (N0*lp_vals_ele)';

end



return
end
    
