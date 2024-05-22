function [x_nodes] = LSC_LP(Bigmat, fvec, numel, nnpe, conn, ...
    totnod, x)
% LSC_angdist smooths the angles after finding the optimum configuration
% from the LP library






x_nodes = zeros(totnod,4);

for i = 1:1:4

    F = zeros(totnod,1);
    

    for iele = 1:1:numel

        fe = zeros(nnpe,1);

        fe(1:nnpe,1) = fvec(iele,1:nnpe)*x(4*iele-4+i);

        [F] = LSC_AssembleFvec(iele,nnpe,conn,fe,F);

    end

    x_nodes(:,i) = Bigmat \ F;



end

