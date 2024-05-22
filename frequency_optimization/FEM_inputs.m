% This script generates the FEA mesh and defines boundary conditions for a
% simply-supported plate with an out-of-plane distributed load appplied in
% negative z-direction (case 1)

% Mesh Inputs



% Total thickness of the plate (assumed to be constant!)
t = 0.01; % meters

Ly = 1;
ratio = 1.5;

nx = 15; 
ny = ceil(ratio*nx);
Lx = Ly/ratio; 

sx = Lx/nx; 
sy = Ly/ny; % size of each element in x and y direction
%%%%%%%%%%%%%%%%%%%%%%% Assign node coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% hold on
crds = zeros((nx+1)*(ny+1),2);
dum = 0;
for j = 0:1:ny
    for i = 0:1:nx
        dum = dum+1;
        crds(dum,1:2) = [i*sx,j*sy];
%         plot(crds(dum,1),crds(dum,2),'*')
    end
end
%%%%%%%%%%%%%%%%%%%%%% Find element connectivity %%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% hold on
% axis equal
conn = zeros(nx*ny, 4);
dum = 0;
for j = 1:1:ny
    for i = 1:1:nx
        dum = dum+1;
        conn(dum,1:1:4) = [i+(j-1)*(nx+1), i+1+(j-1)*(nx+1), ...
            i+1+j*(nx+1), i+j*(nx+1)];
    end
end  


% nnpe: nodes per element
nnpe = 4;

% Number of elements in the mesh
numel = size(conn,1);
% Total number of nodes in the mesh
numnod = size(crds,1);





% BENDING STIFFNESS

% dof: total number of dof per node
dof = 4;

% numind: total number of independent degrees of freedom
numind = max(max(conn))*dof;


% Eltype for bending stiffness
eltyp='BFS';

% Initialize variables
% stiffness matrix
[K, F] =  FEM_initialize(numind);
f = zeros(nnpe*dof,1);

% mass matrix
[M] =  FEM_initialize(numind);






% Force Inputs (point forces and moments on nodes)
% p = 0;
% Fpl_NodComp = [0 0]; % No load in y direction on the upper edge
% [f] = FEM_fpl(f, dof, Fpl_NodComp, p); % Assemble the load vector

% 
% 
% qL_u = -1  ; % in-plane force
% Fipedlu_ElNod = zeros(numel,nnpe);
% 
% % Apply force on the right edge
% Fipedlu_ElNod(nx:nx:nx*nx,2) = 1;
% 
% 
% 
% 
% 
% qL_v = 0; % in-plane force
% Fipedlv_ElNod = zeros(numel,nnpe);


% 
% % BC Inputs (simply-supported ends) 
% %
% % (bottom-left corner fixed, others with w=0 )
% % For the corner nodes 
% dofBCs = [1, 2, 3, 4, ...
%     dof*(nx+1)-dof+1, dof*(nx+1)-dof+2, dof*(nx+1)-dof+3, dof*(nx+1)-dof+4,...
%     dof*(nx+1)*(ny+1)-dof+1, dof*(nx+1)*(ny+1)-dof+2, dof*(nx+1)*(ny+1)-dof+3,  dof*(nx+1)*(ny+1)-dof+4,...
%     dof*(ny*(nx+1)+1)-dof+1, dof*(ny*(nx+1)+1)-dof+2, dof*(ny*(nx+1)+1)-dof+3,  dof*(ny*(nx+1)+1)-dof+4];
% 
% 
% % Along bottom edge: simply support (w=0)
% for i = 2:1:nx
%     dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2];
% end
% 
% % Along left edge: simply support (w=0)
% for i = nx+2:nx+1:(ny-1)*(nx+1)+1
%     dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+3];
% end
% 
% % Along top edge:  simply support (w=0)
% for i = (nx+1)*ny+2:1:(nx+1)*(ny+1)-1
%     dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2];
% end
% 
% % Along right edge: simply support (w=0)
% for i = 2*nx+2:nx+1:(nx+1)*ny
%     dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+3];
% end
% 
% numdofs = size(dofBCs,2); % The total number of BCs
% 
% valdofs = zeros(1,numdofs); % The values of the specified BCs
% 



% BC Inputs (fully-clamped ends) 
%
% (bottom-left corner fixed, others with w=0 )
% For the corner nodes 
dofBCs = [1, 2, 3, 4, ...
    dof*(nx+1)-dof+1, dof*(nx+1)-dof+2, dof*(nx+1)-dof+3, dof*(nx+1)-dof+4,...
    dof*(nx+1)*(ny+1)-dof+1, dof*(nx+1)*(ny+1)-dof+2, dof*(nx+1)*(ny+1)-dof+3,  dof*(nx+1)*(ny+1)-dof+4,...
    dof*(ny*(nx+1)+1)-dof+1, dof*(ny*(nx+1)+1)-dof+2, dof*(ny*(nx+1)+1)-dof+3,  dof*(ny*(nx+1)+1)-dof+4];


% Along bottom edge: simply support (w=0)
for i = 2:1:nx
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3, dof*(i-1)+4];
end

% Along left edge: simply support (w=0)
for i = nx+2:nx+1:(ny-1)*(nx+1)+1
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3, dof*(i-1)+4];
end

% Along top edge:  simply support (w=0)
for i = (nx+1)*ny+2:1:(nx+1)*(ny+1)-1
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3, dof*(i-1)+4];
end

% Along right edge: simply support (w=0)
for i = 2*nx+2:nx+1:(nx+1)*ny
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3, dof*(i-1)+4];
end

numdofs = size(dofBCs,2); % The total number of BCs

valdofs = zeros(1,numdofs); % The values of the specified BCs










% Plot scale
scale=1;

%clear L_x ny ratio nx sx sy
