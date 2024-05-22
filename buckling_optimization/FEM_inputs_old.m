% This script generates the FEA mesh and defines boundary conditions for a
% simply-supported plate with an out-of-plane distributed load appplied in
% negative z-direction (case 1)

% Mesh Inputs



% Total thickness of the plate (assumed to be constant!)
t = 0.127*8; % mm
% t = 8*0.16;

L_x = 254;
% L_x=400;
ny = 16; 
ratio = 1; % x/y
L_y = L_x/ratio;
nx = ratio*ny; 
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
dof_b = 4;

% numind: total number of independent degrees of freedom
numind_b = max(max(conn))*dof_b;


% Eltype for bending stiffness
eltyp_b='BFS';

% Initialize variables
[K_b, F_b] =  FEM_initialize(numind_b);
f_b = zeros(nnpe*dof_b,1);






% MEMBRANE STIFFNESS

% dof: total number of dof per node
dof_m = 2;

% numind: total number of independent degrees of freedom
numind_m = max(max(conn))*dof_m;




% Eltype for bending stiffness
eltyp_m='REC';

% Initialize variables
[K_m, F_m] =  FEM_initialize(numind_m);
f_m = zeros(nnpe*dof_m,1);

% GEOMETRIC STIFFNESS

% Initialize variables
[K_g, F_g] =  FEM_initialize(numind_b);
f_g = zeros(nnpe*dof_b,1);






% Force Inputs (point forces and moments on nodes)
% p = 0;
% Fpl_NodComp = [0 0]; % No load in y direction on the upper edge
% [f] = FEM_fpl(f, dof, Fpl_NodComp, p); % Assemble the load vector



qL_u = 0  ; % in-plane force
Fipedlu_ElNod = zeros(numel,nnpe);

% % Apply force on the right edge
% % Fipedlu_ElNod(nx:nx:nx*nx,2) = 1;
% Fipedlu_ElNod(nx*(ny-1)+1:1:nx*ny,3) = 1;




qL_v = -1; % in-plane force
Fipedlv_ElNod = zeros(numel,nnpe);

% Apply force on the top edge
Fipedlv_ElNod(nx*(ny-1)+1:1:nx*ny,3) = 1;




% BC Inputs (simply-supported ends) 
%

dofBCs_b = [];
% (bottom-left corner fixed, others with w=0 )
% % For the corner nodes 
% dofBCs_b = [1, 2, 3, 4, ...
%     dof_b*(nx+1)-dof_b+1, dof_b*(nx+1)-dof_b+2, dof_b*(nx+1)-dof_b+3, dof_b*(nx+1)-dof_b+4,...
%     dof_b*(nx+1)*(ny+1)-dof_b+1, dof_b*(nx+1)*(ny+1)-dof_b+2, dof_b*(nx+1)*(ny+1)-dof_b+3,  dof_b*(nx+1)*(ny+1)-dof_b+4,...
%     dof_b*(ny*(nx+1)+1)-dof_b+1, dof_b*(ny*(nx+1)+1)-dof_b+2, dof_b*(ny*(nx+1)+1)-dof_b+3,  dof_b*(ny*(nx+1)+1)-dof_b+4];


% Along bottom edge: fixed (w=0)
for i = 1:1:nx+1
    dofBCs_b = [dofBCs_b, dof_b*(i-1)+1, dof_b*(i-1)+2, dof_b*(i-1)+3, dof_b*(i-1)+4 ];
end

% Along left edge: simply support (w=0)
for i = nx+2:nx+1:(ny-1)*(nx+1)+1
    dofBCs_b = [dofBCs_b, dof_b*(i-1)+1, dof_b*(i-1)+3];
end

% Along top edge:  fixed (w=0)
for i = (nx+1)*ny+1:1:(nx+1)*(ny+1)
    dofBCs_b = [dofBCs_b, dof_b*(i-1)+1, dof_b*(i-1)+2, dof_b*(i-1)+3, dof_b*(i-1)+4 ];
end

% Along right edge: simply support (w=0)
for i = 2*nx+2:nx+1:(nx+1)*ny
    dofBCs_b = [dofBCs_b, dof_b*(i-1)+1, dof_b*(i-1)+3];
end

numdofs_b = size(dofBCs_b,2); % The total number of BCs

valdofs_b = zeros(1,numdofs_b); % The values of the specified BCs






% Membrane BCs

dofBCs_m=[];
% % Along left edge: fixed - x (x=0)
% for i = nx+2:nx+1:(ny-1)*(nx+1)+1
%     dofBCs_m = [dofBCs_m, dof_m*(i-1)+1];
% end
% 
% % Along right edge: fixed - x (x=0)
% for i = 2*nx+2:nx+1:ny*(nx+1)
%     dofBCs_m = [dofBCs_m, dof_m*(i-1)+1];
% end

% Along bottom edge: fixed - x & y (y=0)
dofBCs_m_bottom=[];
for i = 1:1:nx+1
    dofBCs_m = [dofBCs_m, dof_m*(i-1)+1, dof_m*(i-1)+2];
    dofBCs_m_bottom = [dofBCs_m_bottom, dof_m*(i-1)+2];
end


% Along top edge:  fixed - x (y=L)
dofBCs_m_top = [];
for i = (nx+1)*ny+1:1:(nx+1)*(ny+1)
    dofBCs_m = [dofBCs_m, dof_m*(i-1)+1];
    dofBCs_m_top = [dofBCs_m_top, dof_m*(i-1)+2];
end


% % Fix the lower left corner node along x and y-directions (x=y=0)
% dofBCs_m = [dofBCs_m, 1, 2];
% 
% % Fix the upper left corner node along x direction (x=0)
% dofBCs_m = [dofBCs_m,  dof_m*(ny*(nx+1)+1)-dof_m+1];



% % Along left edge: fixed - x (x=0)
% for i = 1:nx+1:ny*(nx+1)+1
%     dofBCs_m = [dofBCs_m, dof_m*(i-1)+1];
% end
% 
% % Along bottom edge: fixed - y (y=0)
% for i = 1:1:nx+1
%     dofBCs_m = [dofBCs_m, dof_m*(i-1)+2];
% end




numdofs_m=size(dofBCs_m,2);
valdofs_m=zeros(1,numdofs_m);


% % Top vertical Displacement
% % Along top edge:  fixed - y (y=L)
% for i = (nx+1)*ny+1:1:(nx+1)*(ny+1)
%     dofBCs_m = [dofBCs_m, dof_m*(i-1)+2];
%     valdofs_m = [valdofs_m, -1];
% end




numdofs_m=size(dofBCs_m,2);




% Plot scale
scale=1;

%clear L_x ny ratio nx sx sy
