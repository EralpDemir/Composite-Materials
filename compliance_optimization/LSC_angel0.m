function  [angPly_c] = LSC_angel0(numel, nnpe, conn, nPly, angPlyNodes, N0)
%LSC_angel0 calculate the angles at the element centers

angPly_c = zeros(numel,nPly);

for iPly = 1:1:nPly
    
    ang_vals(:,1) = angPlyNodes(:,iPly);
    
    for iele = 1:1:numel

        ang_vals_ele = zeros(nnpe,1);
        
        for inod = 1:1:nnpe

            NodeNo = conn(iele,inod);

            ang_vals_ele(inod,1) = ang_vals(NodeNo);

        end

        angPly_c(iele,iPly) = N0*ang_vals_ele;

    end
        
end
    
