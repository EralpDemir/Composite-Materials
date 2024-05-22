% This script generates the FEA mesh and defines boundary conditions for a
% plate clamped in all edges with an out-of-plane distributed load appplied
% in negative z-direction (case 3)

% Mesh Inputs

% Eltype
eltyp='BFS';

% Total thickness of the plate (assumed to be constant!)
t = 0.1; % meter

L_x = 1;
ny = 10; 
ratio = 1;
nx = ny*ratio; 
sx = L_x/nx; sy = L_x/ratio/ny; % size of each element in x and y direction
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

conn = zeros(nx*ny, 4);
dum = 0;
for j = 1:1:ny
    for i = 1:1:nx
        dum = dum+1;
        conn(dum,1:1:4) = [i+(j-1)*(nx+1), i+1+(j-1)*(nx+1), ...
            i+1+j*(nx+1), i+j*(nx+1)];
    end
end  
dof = 4; % 1 translation, 2 rotations, 1 twist
nnpe = size(conn,2);
numind = max(max(conn))*dof;
numel = size(conn,1);

% Total number of nodes in the mesh
numnod = size(crds,1);

% Specifications dependent on the element type
% Number quadrature points per element
nqpts=0;

% Initialize variables
[K, f, totdof, dofst, dofID, ndir, numenod, numvnod ] = ...
    FEM_initialize(conn,numel,numnod,eltyp,numind);

% Force Inputs (point forces and moments on nodes)
p = 0;
Fpl_NodComp = [0 0]; % No load in y direction on the upper edge
[f] = FEM_fpl(f, dof, Fpl_NodComp, p); % Assemble the load vector

qS_w = -1; % out-of-plane force
Fsdlw_El = ones(numel,1);

qL_w = 0; % out-of-plane force
Fedlw_ElNod = zeros(numel,nnpe);

qL_u = 0; % in-plane force
Fipedlu_ElNod = zeros(numel,nnpe);

qL_v = 0; % in-plane force
Fipedlv_ElNod = zeros(numel,nnpe);

% BC Inputs (fixed ends - all DOFs are zero on all edges)
%
% For the corner nodes:
dofBCs = [1, 2, 3, 4, ...
        dof*(nx+1)-3, dof*(nx+1)-2, dof*(nx+1)-1, dof*(nx+1), ...
        dof*(nx*(ny+1)+1)-3, dof*(nx*(ny+1)+1)-2, dof*(nx*(ny+1)+1)-1, dof*(nx*(ny+1)+1), ...
        dof*(nx+1)*(ny+1)-3, dof*(nx+1)*(ny+1)-2, dof*(nx+1)*(ny+1)-1, dof*(nx+1)*(ny+1)];
    
% Along bottom edge: Fixed support (u=0,v=0,w=0,theta_x=0,theta_y=0,wxy=0)
for i = 2:1:nx
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3 ];
end

% Along left edge: Fixed support (u=0,v=0,w=0,theta_x=0,theta_y=0,wxy=0)
for i = nx+2:nx+1:(ny-1)*(nx+1)+1
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3];
end

% Along top edge: Fixed support (u=0,v=0,w=0,theta_x=0,theta_y=0,wxy=0)
for i = (nx+1)*ny+2:1:(nx+1)*(ny+1)-1
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3];
end

% Along right edge: Fixed support (u=0,v=0,w=0,theta_x=0,theta_y=0,wxy=0)
for i = 2*nx+2:nx+1:(nx+1)*ny
    dofBCs = [dofBCs, dof*(i-1)+1, dof*(i-1)+2, dof*(i-1)+3];
end

numBCs = size(dofBCs,2); % The total number of BCs

valBCs = zeros(1,numBCs); % The values of the specified BCs

% Plot scale
scale=100;

