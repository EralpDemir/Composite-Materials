function [K, f, totdof, dofst, dofID, ndir, numenod, numvnod ] = ...
    FEM_initialize(conn,numel,numnod,eltyp,numind)
%FEM_initialize takes the element type and mesh nodal information and
%initialized the FEA matrices related to each element type

switch eltyp
    
    case 'TRI' %__________________________________________________________
        
        
        ndir = 0;
        dofst = 0;
        dofID = 0;
        numenod = 0;
        numvnod = 0;
        totdof = numind;
        
        K = zeros(numind);
        f = zeros(numind,1);
        
    case 'HCT' %__________________________________________________________
        
        %%%%%%%%%%%%%%%%%%%%%%%% mark the global vertex and edge nodes
        for i = 1:numel
            dofID(conn(i,1),:) = 0; dofID(conn(i,2),:) = 0; dofID(conn(i,3),:) = 0;  % vertex nodes
            dofID(conn(i,4),:) = 1; dofID(conn(i,5),:) = 1; dofID(conn(i,6),:) = 1;  % edge nodes
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL NODES FREEDOM
        % Starting dof indicator for each node respect to global node sequentially
        for i = 1:numnod
            dofst(i,:) = 1;
            for j=1:i-1
                if(dofID(j) == 1)    % edge node,   1 d.o.f.
                    dofst(i,:) = dofst(i,:)+1;
                else                % vertex node, 3 d.o.f.
                    dofst(i,:) = dofst(i,:)+3;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% COUNT SYSTEM TOTAL GLOBAL NODE
        % find the coordinates of the edge nodes with respect to global nodes
        numenod = 0;  % number of edge nodes
        for i = 1:numnod  % number of global nodes
            if(dofID(i) == 1)
                numenod = numenod+1;
            end
        end
        
        numvnod = 0 ; % number of vertex nodes
        for i = 1:numnod
            if(dofID(i) == 0)
                numvnod = numvnod+1;
            end
        end
        
        % total number of global modes
        totdof = 3*numvnod+numenod;  % number of global modes
        
        % Stiffness matrix
        K = zeros(totdof);
        f = zeros(totdof,1);
        
        % initialize the orientation index of the normal derivative at the midside
        % nodes ndir is 9*1
        for i = 1:numnod
            ndir(i,:) = 1.0;
        end
        
    case 'REC' %__________________________________________________________
        
        ndir = 0;
        dofst = 0;
        dofID = 0;
        K = zeros(numind);
        f = zeros(numind,1);
        numenod = 0;
        numvnod = 0;
        totdof = numind;
        
    case 'BFS' %__________________________________________________________
        
        ndir = 0;
        dofst = 0;
        dofID = 0;
        K = zeros(numind);
        f = zeros(numind,1);
        numenod = 0;
        numvnod = 0;
        totdof = numind;
        
end

