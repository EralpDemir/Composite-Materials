function [lpNodes] = LSC_lpdist(BigNmat, Bigfvec, numel, nnpe, conn, ...
    totnod,xi)
% LSC_angdist smooths the angles after finding the optimum configuration
% from the LP library

lpNodes = zeros(totnod,1);



F = zeros(totnod,1);

for iele = 1:1:numel

    fe = zeros(nnpe,1);

    fe(1:nnpe,1) = Bigfvec(iele,1:nnpe)*xi(iele);

    [F] = LSC_AssembleFvec(iele,nnpe,conn,fe,F);

end

lpNodes = BigNmat \ F;




return

end

