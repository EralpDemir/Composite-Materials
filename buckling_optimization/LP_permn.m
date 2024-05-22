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
if nPly==1
    N = 1;
else
    N = nPly/2;
end
l = length(ang_range);
result = zeros(l^N,N);
for i = 1:N
    a = reshape(reshape(repmat((1:l),1,l^(N-1)),l^i,l^(N-i))',l^N,1);
    result(:,i) = ang_range(a);
end


