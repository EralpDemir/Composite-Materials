% Eralp Demir
% 02.12.2019
% Updated 28.02.2021

clear
close all
clc

% Load data
load data



% Load mesh data
load mesh
% 
% % load BCs
% load BCs


% Inputs of LSC
inputs


% Initialize LSC
initialize


% Plot the angles
outputel(conn,crds,numel,nnpe,All_angPly_c(:,izPly))



% get the trend (0 or 90 for plates)
av_ang = mean(All_angPly_c(:,izPly));
if av_ang<=45 && av_ang>=-45
    av_ang=0;
else
    av_ang=90;
end


% av_ang=0;


% Based on the trend inflow BCs

% If horizontal flow ==> assign left & right surfaces as zero flow BC

if av_ang==0
    % Left end
    val_x = min(crds(:,1));
    nodes = find(crds(:,1)==val_x);
    dofBCs = nodes;
    numBCs = size(nodes,1);

    
    % Right end
    val_x = max(crds(:,1));
    nodes = find(crds(:,1)==val_x);
    dofBCs = [dofBCs; nodes];
    numBCs = numBCs + size(nodes,1);
    
    valBCs = zeros(1,numBCs);
    
elseif av_ang==90
    
    % Bottom end
    val_y = min(crds(:,2));
    nodes = find(crds(:,2)==val_y);
    dofBCs = nodes;
    numBCs = size(nodes,1);

    
    % Top end
    val_y = max(crds(:,2));
    nodes = find(crds(:,2)==val_y);
    dofBCs = [dofBCs; nodes];
    numBCs = numBCs + size(nodes,1);
    
    valBCs = zeros(1,numBCs);
    
    
end

% If vertical flow ==> assign left & right surfaces as zero flow BC


% 
% % Data for Tom's case
% angles=angles*pi/180;

% % WITH LSC
% angles(1:1:totnod,1) = angPlyNodes(:,izPly)*pi/180;


% WITHOUT LSC
% Process angles to find the angles at the nodes
angles=zeros(totnod,1);
for inod=1:1:totnod

    counter=0;
    for j=1:1:numel
        
        for k=1:1:nnpe
        
            if conn(j,k)==inod
                counter = counter + 1;
                angles(inod)=angles(inod)+All_angPly_c(j,izPly); 
            end
        end
    end
    
    angles(inod)=angles(inod)/counter*pi/180;
    
    
                
        
    
end

% Plot nodal angles
outputnod(crds,conn,nnpe,numel,angles)


% 
% 
% % angles=All_angPly(:,izPly);




% % Calculate BigNmat
% calc_BigMmat

% Calculate BigNmat
calc_BigNmat

% Overall stiffness matrix
BigKmat = BigNmat + alpha * BigCmat + beta * BigDmat;


% Apply boundary conditions
[BigKmat, Bigfvec] = applybcs(BigKmat, Bigfvec, dofBCs, valBCs, numBCs);
% Only at node=1, thickness=1
% BigKmat=BigNmat + BigMmat;
% dofBCs=121;
% numBCs=size(dofBCs,1);
% valBCs=zeros(numBCs,1);

% rank(BigNmat)

% [BigNmat, Bigfvec] = applybcs(BigNmat, Bigfvec, dofBCs, valBCs, numBCs);

% rank(BigNmat)

% Find the thickness distribution
thick = -BigKmat\Bigfvec;




% thick=thick/max(abs(thick))*exp(1);


% % Plot thickness (nodal)
% output_thickness(crds,totnod,thick);


% Find elemental thickness (elemental)
[thick_el] = thick_el_center(numel, nnpe, conn, thick, N0);





% Stream function calculation
analyze







