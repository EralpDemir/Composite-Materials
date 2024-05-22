function [angPlyNodes] = LSC_angdist(Bigmat, fvec, numel, nnpe, conn, ...
    totnod, nPly, angPly)
% LSC_angdist smooths the angles after finding the optimum configuration
% from the LP library

angPlyNodes = zeros(totnod,nPly);

for iPly = 1:1:nPly

    F = zeros(totnod,1);

    for iele = 1:1:numel

        fe = zeros(nnpe,1);

        fe(1:nnpe,1) = fvec(iele,1:nnpe)*angPly(iele,iPly);

        [F] = LSC_AssembleFvec(iele,nnpe,conn,fe,F);

    end

    angs = Bigmat \ F;

    angPlyNodes(:,iPly) = angs;

end

