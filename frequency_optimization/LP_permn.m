function result = LP_permn(ang_range, nPly)

%LP_Permn   calculates the lamination parameters space
%
% INPUTS:
% ang_range: range of possible fiber angles in each lamina
% nPly: total number of plies
%
% OUTPUTS:
% result: range of possible angles in each ply
%
N = nPly/2;
l = length(ang_range);
a = zeros(l^N,N);
for i = 1:N
    a(:,i) = reshape(reshape(repmat((1:l),1,l^(N-1)),l^i,l^(N-i))',l^N,1);
end
result = ang_range(a);

