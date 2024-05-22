function [stress_crds, max_stress_vals, max_strs_els] = ...
    FEM_output_stressdist(nPly,numel,conn,crds,nnpe,dof,dofst,numind,...
    d,Q_bar,angPly,zPly,eltyp)
%FEM_output_stressdist plots the stress distribution and generates the
%location, values, and elements number of maximum stresses

fig_num = 0;
Estrain = zeros(3,numel);
Estress = zeros(3,numel);

switch eltyp
    %_____________________________ REC _____________________________
    case 'REC'
        
        %loop through each ply
        for j = 1:1:nPly
            for ielp = 1:numel
                nd = conn(ielp,1:nnpe); % Loop for total number of elements
                
                [index] = FEM_elementdof(nd,nnpe,dof,eltyp); % extract system dofs associated with element
                
                eldisp = d(index);
                
                % strain-displacement matrix for the 2a*2b  rectangular element
                a = abs(crds(conn(ielp,2),1)-crds(conn(ielp,1),1))/2;
                b = abs(crds(conn(ielp,4),2)-crds(conn(ielp,1),2))/2;
                ksi = 0;   % center of an element
                eta = 0;   % center of an element
                
                % define strain displacement and find the value in center
                B = [-(b-ksi)/4/a/b  0  (b-ksi)/4/a/b  0  (b+ksi)/4/a/b  0  -(b+ksi)/4/a/b  0;...
                    0 -(a-eta)/4/a/b  0  -(a+eta)/4/a/b  0  (a+eta)/4/a/b  0  (a-eta)/4/a/b;...
                    -(a-eta)/4/a/b -(b-ksi)/4/a/b -(a+eta)/4/a/b (b-ksi)/4/a/b ...
                    (a+eta)/4/a/b (b+ksi)/4/a/b (a-eta)/4/a/b -(b+ksi)/4/a/b];
                
                q = Q_bar(:,:,j);               % reduced stiffness matrix for each ply
                Estrain(:,ielp) = B*eldisp ;              % Compute Strains
                Estress(:,ielp) = q*B*eldisp ;             % Compute Stresses
                
            end
            strain(:,:,j) = Estrain;
            stress(:,:,j) = Estress;
        end
        M = []; max_stress_vals = [];
        for j = 1:3
            Z = [];
            [~,m] = find(abs(stress) == max(max(abs(stress(j,:,:)))));
            ind = mod(max(m),numel);
            if ind == 0
                ind = numel;
            end
            M = [M ind];
            
            EPS = [];
            SIG = [];
            for i = 1:nPly
                z = linspace(zPly(ind,i,1), zPly(ind,i,2), 2);
                Z = [Z z];
                EPS = [EPS strain(j,ind,i)*ones(1,2)];
                SIG = [SIG stress(j,ind,i)*ones(1,2)];
            end
            
            Mask = [];
            for  i = (nPly/2)*2+1 : 2 : nPly*2
                Mask = [ i:i+2-1 Mask];
            end
            Mask = [Mask 1:nPly*2/2];
            
            %%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(EPS(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((EPS(i)+EPS(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Strain distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\epsilon_',num2str(j),' [mm]'])
            ylabel('laminate thickness [m]')
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(SIG(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((SIG(i)+SIG(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Stress distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\sigma_',num2str(j),' [MPa]'])
            ylabel('laminate thickness [m]')
            
            stress_crds(:,j) = [sum(crds(conn(ind,:),1))/size(conn,2); ...
                sum(crds(conn(ind,:),2))/size(conn,2)];
            AA = sort(SIG(Mask));
            if  abs(AA(1)) < abs(AA(end))
                MSV = max(abs(SIG(Mask)));
            else
                MSV = -max(abs(SIG(Mask)));
            end
            max_stress_vals = [max_stress_vals MSV];
        end
        max_strs_els = M;
        
        %_____________________________ BFS _____________________________
    case 'BFS'
        
        Estrain = zeros(3,numel);
        Estress = zeros(3,numel);
        Z = [];
        Strain = []; Stress = [];
        %loop through each ply
        for j = 1:1:nPly
            z = linspace(zPly(1,j,1), zPly(1,j,2), 2);
            Z = [Z z];
            for k = 1:length(z)
                for ielp = 1:numel
                    nd = conn(ielp,1:nnpe); % Loop for total number of elements
                    
                    [index] = FEM_elementdof(nd,nnpe,dof,eltyp); % extract system dofs associated with element
                    
                    eldisp = d(index);
                    
                    % strain-displacement matrix for the 2a*2b  BFS element
                    a = abs(crds(conn(ielp,2),1)-crds(conn(ielp,1),1))/2;
                    b = abs(crds(conn(ielp,4),2)-crds(conn(ielp,1),2))/2;
                    ksi = 0;   % center of an element
                    eta = 0;   % center of an element
                    
                    % define strain displacement and find the value in center
                    B = [                      ((3*ksi*eta^3)/8 - (9*ksi*eta)/8 + (3*ksi)/4)/a^2,                   ((3*eta)/8 + (3*ksi)/4 - (9*eta*ksi)/8 + (3*eta^3*ksi)/8 - eta^3/8 - 1/4)/a,                       (b*((3*ksi*eta^3)/8 - (3*ksi*eta^2)/8 - (3*ksi*eta)/8 + (3*ksi)/8))/a^2,               (b*(eta/8 + (3*ksi)/8 - (3*eta*ksi)/8 - (3*eta^2*ksi)/8 + (3*eta^3*ksi)/8 + eta^2/8 - eta^3/8 - 1/8))/a,                    -((3*ksi*eta^3)/8 - (9*ksi*eta)/8 + (3*ksi)/4)/a^2,                      ((3*ksi)/4 - (3*eta)/8 - (9*eta*ksi)/8 + (3*eta^3*ksi)/8 + eta^3/8 + 1/4)/a,                       -(b*((3*ksi*eta^3)/8 - (3*ksi*eta^2)/8 - (3*ksi*eta)/8 + (3*ksi)/8))/a^2,              -(b*(eta/8 - (3*ksi)/8 + (3*eta*ksi)/8 + (3*eta^2*ksi)/8 - (3*eta^3*ksi)/8 + eta^2/8 - eta^3/8 - 1/8))/a,                     -((3*ksi)/4 + (9*eta*ksi)/8 - (3*eta^3*ksi)/8)/a^2,                     ((3*eta)/8 + (3*ksi)/4 + (9*eta*ksi)/8 - (3*eta^3*ksi)/8 - eta^3/8 + 1/4)/a,                         (b*((3*ksi)/8 + (3*eta*ksi)/8 - (3*eta^2*ksi)/8 - (3*eta^3*ksi)/8))/a^2,              -(b*(eta/8 + (3*ksi)/8 + (3*eta*ksi)/8 - (3*eta^2*ksi)/8 - (3*eta^3*ksi)/8 - eta^2/8 - eta^3/8 + 1/8))/a,                     ((3*ksi)/4 + (9*eta*ksi)/8 - (3*eta^3*ksi)/8)/a^2,                   -((3*eta)/8 - (3*ksi)/4 - (9*eta*ksi)/8 + (3*eta^3*ksi)/8 - eta^3/8 + 1/4)/a,                         -(b*((3*ksi)/8 + (3*eta*ksi)/8 - (3*eta^2*ksi)/8 - (3*eta^3*ksi)/8))/a^2,               (b*(eta/8 - (3*ksi)/8 - (3*eta*ksi)/8 + (3*eta^2*ksi)/8 + (3*eta^3*ksi)/8 - eta^2/8 - eta^3/8 + 1/8))/a;...
                        ((3*eta*ksi^3)/8 - (9*eta*ksi)/8 + (3*eta)/4)/b^2,                       (a*((3*eta*ksi^3)/8 - (3*eta*ksi^2)/8 - (3*eta*ksi)/8 + (3*eta)/8))/b^2,                   ((3*eta)/4 + (3*ksi)/8 - (9*eta*ksi)/8 + (3*eta*ksi^3)/8 - ksi^3/8 - 1/4)/b,               (a*((3*eta)/8 + ksi/8 - (3*eta*ksi)/8 - (3*eta*ksi^2)/8 + (3*eta*ksi^3)/8 + ksi^2/8 - ksi^3/8 - 1/8))/b,                     ((3*eta)/4 + (9*eta*ksi)/8 - (3*eta*ksi^3)/8)/b^2,                         -(a*((3*eta)/8 + (3*eta*ksi)/8 - (3*eta*ksi^2)/8 - (3*eta*ksi^3)/8))/b^2,                    ((3*eta)/4 - (3*ksi)/8 + (9*eta*ksi)/8 - (3*eta*ksi^3)/8 + ksi^3/8 - 1/4)/b,              -(a*((3*eta)/8 - ksi/8 + (3*eta*ksi)/8 - (3*eta*ksi^2)/8 - (3*eta*ksi^3)/8 + ksi^2/8 + ksi^3/8 - 1/8))/b,                     -((3*eta)/4 + (9*eta*ksi)/8 - (3*eta*ksi^3)/8)/b^2,                         (a*((3*eta)/8 + (3*eta*ksi)/8 - (3*eta*ksi^2)/8 - (3*eta*ksi^3)/8))/b^2,                     ((3*eta)/4 + (3*ksi)/8 + (9*eta*ksi)/8 - (3*eta*ksi^3)/8 - ksi^3/8 + 1/4)/b,              -(a*((3*eta)/8 + ksi/8 + (3*eta*ksi)/8 - (3*eta*ksi^2)/8 - (3*eta*ksi^3)/8 - ksi^2/8 - ksi^3/8 + 1/8))/b,                    -((3*eta*ksi^3)/8 - (9*eta*ksi)/8 + (3*eta)/4)/b^2,                       -(a*((3*eta*ksi^3)/8 - (3*eta*ksi^2)/8 - (3*eta*ksi)/8 + (3*eta)/8))/b^2,                      ((3*eta)/4 - (3*ksi)/8 - (9*eta*ksi)/8 + (3*eta*ksi^3)/8 + ksi^3/8 + 1/4)/b,               (a*((3*eta)/8 - ksi/8 - (3*eta*ksi)/8 - (3*eta*ksi^2)/8 + (3*eta*ksi^3)/8 - ksi^2/8 + ksi^3/8 + 1/8))/b;...
                        -(2*(- (9*eta^2*ksi^2)/16 + (9*eta^2)/16 + (9*ksi^2)/16 - 9/16))/(a*b), (2*((9*eta^2*ksi^2)/16 - (3*eta^2*ksi)/8 - (3*eta^2)/16 - (9*ksi^2)/16 + (3*ksi)/8 + 3/16))/b, (2*((9*eta^2*ksi^2)/16 - (9*eta^2)/16 - (3*eta*ksi^2)/8 + (3*eta)/8 - (3*ksi^2)/16 + 3/16))/a, (9*eta^2*ksi^2)/8 - (3*eta^2*ksi)/4 - (3*eta^2)/8 - (3*eta*ksi^2)/4 + (eta*ksi)/2 + eta/4 - (3*ksi^2)/8 + ksi/4 + 1/8, (2*(- (9*eta^2*ksi^2)/16 + (9*eta^2)/16 + (9*ksi^2)/16 - 9/16))/(a*b), -(2*(- (9*eta^2*ksi^2)/16 - (3*eta^2*ksi)/8 + (3*eta^2)/16 + (9*ksi^2)/16 + (3*ksi)/8 - 3/16))/b, -(2*((9*eta^2*ksi^2)/16 - (9*eta^2)/16 - (3*eta*ksi^2)/8 + (3*eta)/8 - (3*ksi^2)/16 + 3/16))/a, (9*eta^2*ksi^2)/8 + (3*eta^2*ksi)/4 - (3*eta^2)/8 - (3*eta*ksi^2)/4 - (eta*ksi)/2 + eta/4 - (3*ksi^2)/8 - ksi/4 + 1/8, -(2*(- (9*eta^2*ksi^2)/16 + (9*eta^2)/16 + (9*ksi^2)/16 - 9/16))/(a*b), (2*(- (9*eta^2*ksi^2)/16 - (3*eta^2*ksi)/8 + (3*eta^2)/16 + (9*ksi^2)/16 + (3*ksi)/8 - 3/16))/b, (2*(- (9*eta^2*ksi^2)/16 + (9*eta^2)/16 - (3*eta*ksi^2)/8 + (3*eta)/8 + (3*ksi^2)/16 - 3/16))/a, (9*eta^2*ksi^2)/8 + (3*eta^2*ksi)/4 - (3*eta^2)/8 + (3*eta*ksi^2)/4 + (eta*ksi)/2 - eta/4 - (3*ksi^2)/8 - ksi/4 + 1/8, (2*(- (9*eta^2*ksi^2)/16 + (9*eta^2)/16 + (9*ksi^2)/16 - 9/16))/(a*b), -(2*((9*eta^2*ksi^2)/16 - (3*eta^2*ksi)/8 - (3*eta^2)/16 - (9*ksi^2)/16 + (3*ksi)/8 + 3/16))/b, -(2*(- (9*eta^2*ksi^2)/16 + (9*eta^2)/16 - (3*eta*ksi^2)/8 + (3*eta)/8 + (3*ksi^2)/16 - 3/16))/a, (9*eta^2*ksi^2)/8 - (3*eta^2*ksi)/4 - (3*eta^2)/8 + (3*eta*ksi^2)/4 - (eta*ksi)/2 - eta/4 - (3*ksi^2)/8 + ksi/4 + 1/8];
                    
                    q = Q_bar(:,:,j);               % reduced stiffness matrix for each ply
                    Estrain(:,ielp) = z(k)*B*eldisp ;              % Compute curvatures
                    Estress(:,ielp) = q*z(k)*B*eldisp ;             % Compute Stresses
                end
                strain(:,:,k) = Estrain;
                stress(:,:,k) = Estress;
            end
            Strain = [Strain strain];
            Stress = [Stress stress];
        end
        
        % find maximum stress in a ply and stress for corresponding element in
        % other plies
        M = []; max_stress_vals = [];
        for j = 1:3
            [~,m] = find(abs(Stress) == max(max(abs(Stress(j,:,:)))));  %index of an element with maximum stress
            ind = mod(mod(max(m),numel*nPly),numel); % find corresponding element index for other layers
            %if index is given zero consider it as numel (we do not have zero index)
            if ind == 0
                ind = numel;
            end
            M = [M ind];
            
            EPS = []; SIG = [];
            for i=1:nPly
                Eps = zeros(1,length(2));
                for k  = 1:2
                    Eps(k) = Strain(j,ind+(i-1)*numel,k);
                    Sig(k) = Stress(j,ind+(i-1)*numel,k);
                end
                EPS = [EPS Eps];
                SIG = [SIG Sig];
            end
            
            Mask = [];
            for  i = (nPly/2)*2+1 : 2 : nPly*2
                Mask = [ i:i+2-1 Mask];
            end
            Mask = [Mask 1:nPly*2/2];
            
            %%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(EPS(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((EPS(i)+EPS(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Strain distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\epsilon_',num2str(j),' [mm]'])
            ylabel('laminate thickness [m]')
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(SIG(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((SIG(i)+SIG(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Stress distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\sigma_',num2str(j),' [MPa]'])
            ylabel('laminate thickness [m]')
            
            stress_crds(:,j) = [sum(crds(conn(ind,:),1))/size(conn,2) ; ...
                sum(crds(conn(ind,:),2))/size(conn,2)];
            AA = sort(SIG(Mask));
            if  abs(AA(1)) < abs(AA(end))
                MSV = max(abs(SIG(Mask)));
            else
                MSV = -max(abs(SIG(Mask)));
            end
            max_stress_vals = [max_stress_vals MSV];
        end
        max_strs_els = M;
        
        %_____________________________ TRI _____________________________
    case 'TRI'
        
        %loop through each ply
        for j = 1:1:nPly
            for ielp = 1:numel
                nd = conn(ielp,1:nnpe); % Loop for total number of elements
                
                [index] = FEM_elementdof(nd,nnpe,dof,eltyp); % extract system dofs associated with element
                
                eldisp = d(index);
                
                x1 = crds(conn(ielp,1),1);  y1 = crds(conn(ielp,1),2);
                
                x2 = crds(conn(ielp,2),1);  y2 = crds(conn(ielp,2),2);
                
                x3 = crds(conn(ielp,3),1);  y3 = crds(conn(ielp,3),2);
                
                b1 = y2-y3;
                c1 = x3-x2;
                
                b2 = y3-y1;
                c2 = x1-x3;
                
                b3 = y1-y2;
                c3 = x2-x1;
                
                area = det([1 x1 y1;
                    1 x2 y2;
                    1 x3 y3])/2;
                
                B = [b1 0 b2 0 b3 0;...
                    0 c1 0 c2 0 c3;...
                    c1 b1 c2 b2 c3 b3]/(2*area);
                
                q = Q_bar(:,:,j);                         % reduced stiffness matrix for each ply
                Estrain(:,ielp) = B*eldisp ;              % Compute Strains
                Estress(:,ielp) = q*B*eldisp ;            % Compute Stresses
                
            end
            strain(:,:,j) = Estrain;
            stress(:,:,j) = Estress;
        end
        
        M = []; max_stress_vals = [];
        for j = 1:3
            Z = [];
            [~,m] = find(abs(stress) == max(max(abs(stress(j,:,:)))));
            ind = mod(max(m),numel);
            if ind == 0
                ind = numel;
            end
            M = [M ind];
            
            EPS = [];
            SIG = [];
            for i = 1:nPly
                z = linspace(zPly(ind,i,1), zPly(ind,i,2), 2);
                Z = [Z z];
                EPS = [EPS strain(j,ind,i)*ones(1,2)];
                SIG = [SIG stress(j,ind,i)*ones(1,2)];
            end
            
            Mask = [];
            for  i = (nPly/2)*2+1 : 2 : nPly*2
                Mask = [ i:i+2-1 Mask];
            end
            Mask = [Mask 1:nPly*2/2];
            
            %%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(EPS(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((EPS(i)+EPS(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Strain distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\epsilon_',num2str(j),' [mm]'])
            ylabel('laminate thickness [m]')
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(SIG(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((SIG(i)+SIG(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Stress distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\sigma_',num2str(j),' [MPa]'])
            ylabel('laminate thickness [m]')
            
            stress_crds(:,j) = [sum(crds(conn(ind,:),1))/size(conn,2) ; sum(crds(conn(ind,:),2))/size(conn,2)];
            AA = sort(SIG(Mask));
            if  abs(AA(1)) < abs(AA(end))
                MSV = max(abs(SIG(Mask)));
            else
                MSV = -max(abs(SIG(Mask)));
            end
            max_stress_vals = [max_stress_vals MSV];
        end
        max_strs_els = M;
        
        %_____________________________ HCT _____________________________
    case 'HCT'
        
        Z = [];
        Strain = []; Stress = [];
        % loop through each ply
        for j = 1:1:nPly
            z = linspace(zPly(1,j,1), zPly(1,j,2), 2);
            Z = [Z z];
            for k = 1:length(z)
                for ielp = 1:numel
                    
                    [index] = FEM_elementdof_hct(conn(ielp,1:nnpe),dofst); % extract system dofs associated with element
                    
                    eldisp = d(index);
                    
                    % Calculate strain-displacement matrix
                    [lpsi] = FEM_ElStrDispB(ielp, crds,conn, numind);
                    
                    q = Q_bar(:,:,j);               % reduced stiffness matrix for each ply
                    Estrain(:,ielp) = z(k)*lpsi*eldisp ;              % Compute curvatures
                    Estress(:,ielp) = q*z(k)*lpsi*eldisp ;             % Compute Stresses
                end
                strain(:,:,k) = Estrain;
                stress(:,:,k) = Estress;
            end
            Strain = [Strain strain];
            Stress = [Stress stress];
            
        end
        
        M = []; max_stress_vals = [];
        for j = 1:3
            [~,m] = find(abs(Stress) == max(max(abs(Stress(j,:,:)))));  %index of an element with maximum stress
            ind = mod(mod(max(m),numel*nPly),numel); % find corresponding element index for other layers
            %if index is given zero consider it as numel (we do not have zero index)
            if ind == 0
                ind = numel;
            end
            M = [M ind];
            
            EPS = []; SIG = [];
            for i = 1:nPly
                Eps = zeros(1,length(2));
                for k  = 1:2
                    Eps(k) = Strain(j,ind+(i-1)*numel,k);
                    Sig(k) = Stress(j,ind+(i-1)*numel,k);
                end
                EPS = [EPS Eps];
                SIG = [SIG Sig];
            end
            
            Mask = [];
            for  i = (nPly/2)*2+1 : 2 : nPly*2
                Mask = [ i:i+2-1 Mask];
            end
            Mask = [Mask 1:nPly*2/2];
            
            %%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(EPS(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((EPS(i)+EPS(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Strain distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\epsilon_',num2str(j),' [mm]'])
            ylabel('laminate thickness [m]')
            
            figure( fig_num + 1 );
            fig_num = fig_num + 1;
            plot(SIG(Mask),Z(Mask))
            for i = 1:2:length(Mask)-1
                nPl = (i+1)/2;
                text((SIG(i)+SIG(i+1))/2, (Z(i)+Z(i+1))/2, ...
                    ['\theta_',num2str(nPl),' = ', num2str(angPly(ind,nPl))])
            end
            %             title( ['Stress distribution along the laminate thickness in direction(',num2str(j),')'] )
            xlabel(['\sigma_',num2str(j),' [MPa]'])
            ylabel('laminate thickness [m]')
            
            stress_crds(:,j) = [sum(crds(conn(ind,1:3),1))/3 ; sum(crds(conn(ind,1:3),2))/3];
            AA = sort(SIG(Mask));
            if  abs(AA(1)) < abs(AA(end))
                MSV = max(abs(SIG(Mask)));
            else
                MSV = -max(abs(SIG(Mask)));
            end
            max_stress_vals = [max_stress_vals MSV];
        end
        max_strs_els = M;
        
end


