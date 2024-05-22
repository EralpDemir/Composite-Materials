% This script defines the thickness, fiber orientations, type of laminate
% and it's specifications

%__________________ Inputs for the layers of the laminate ________________

% nPlys SHALL BE AN EVEN NUMBER FOR SYMMETRIC LAYUP!
nPlys = 8;

%%%%%%%%%%%%%%%%%%%%%%% Define the possible fiber angles
div = 1;
ang_range = -90:div:90-div;
% ang_range = 0:div:180-div;

%%%%%%%%%%%%%%%%%%%%%%% Define the type of laminate
% Single:     1
% Symmetric:  2
% Balanced:   3
% Symmetric & Balanced 4
% Asymmetric: 5
% Specially orthotropic: 7
Ltype = 7;

%%%%%%%%%%%%%%%%%%%%%%% Specification of z-ccodinate of each ply
if Ltype == 1
    nPlys = 1;
end

if nPlys == 1
    zPly(1:numel, 1, 1) = -t/2;
    zPly(1:numel, 1, 2) = t/2;
else
    nPlyV = 0:nPlys-1;
    zPly = zeros(numel, nPlys, 2);
    for i = 1:nPlys/2
        zPly(1:numel, i, 1) = t/nPlys*(nPlyV(i));
        zPly(1:numel, i, 2) = t/nPlys*(nPlyV(i+1));
    end
    zPly(:, :, 1) = [zPly(:, 1:nPlys/2, 1), -zPly(:, 1:nPlys/2, 2)];
    zPly(:, :, 2) = [zPly(:, 1:nPlys/2, 2), -zPly(:, 1:nPlys/2, 1)];
end
clear nPlyV

%--------------------------------------------------------------------------
% The order of layers in the laminates is as following schemes:
%
%       'Even Layer'              'Single Layer'
%    -----------------          -----------------
%           2                   - - - - 1 - - - -
%    -----------------          -----------------
%           1
%    - - - - - - - - -   <==== mid-plane
%           3
%    -----------------
%           4
%    -----------------
%--------------------------------------------------------------------------
