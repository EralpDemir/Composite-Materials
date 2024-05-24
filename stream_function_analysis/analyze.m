

% 
% e = 2; g = 1;
% [x,y] = meshgrid(0:20,0:15);  % This makes regular grid
% u = e*x-g*y;                  % Linear velocity field
% v = g*x-e*y;
% 
% 
% 
% for j=1:21
%      for i=1:16
%          u(i,j)=u(i,j)/sqrt(u(i,j)^2+v(i,j)^2);
%          v(i,j)=v(i,j)/sqrt(u(i,j)^2+v(i,j)^2);
%      end
% end

t = thick;
x = crds(:,1);
y = crds(:,2);

norm_U = exp(t);
u = norm_U.*cos(angles);
v = norm_U.*sin(angles);


x_seeds = unique(x);
y_seeds = unique(y);

[X,Y] = meshgrid(unique(x),unique(y));

U = griddata(x,y,u,X,Y);
V = griddata(x,y,v,X,Y);
th = griddata(x,y,norm_U,X,Y);

[dU_dx,dU_dy] = gradient(U);
[dV_dx,dV_dy] = gradient(V);

% Check continuity
figure;
divU = dU_dx + dV_dy;
contourf(X,Y,divU);
colorbar;
title('div(u)')
axis equal
% % Random angles
% for i=1:20
%     for j=1:15
%         th=pi/2*rand()-pi/2;
%         u(i,j)=cos(th);
%         v(i,j)=sin(th);
%     end
% end

figure;
contourf(X,Y,th);
colorbar;
title('thickness')
axis equal


[psi] = flowfun(U,V,'-');  % Here comes the potential and streamfun.

%  contour(phi,20,'--r')   % Contours of potential
figure
hold on
contour(X,Y,psi, 40)    % Contours of streamfunction
title('stream function')
axis equal

figure
hold on
quiver(X,Y,U,V,'k')         % Now superimpose the velocity field
title('velocity')
axis equal

% % Plot the angles
% scale=5e-3;
% for i=1:1:size(x,1)
%     line([x(i), x(i)+scale*cos(theta(i))], ...
%         [y(i), y(i)+scale*sin(theta(i))],'Color','r');
% end


% sgnU=sign(mean(mean(U)));
% sgnV=sign(mean(mean(V)));



% Plot the streamlines
figure
hold on

% Average angles reveal starting edge
if av_ang==0
    sgnU=sign(mean(mean(U)));
    sgnV=sgnU;
    % Start from left edge
    st_x = min(crds(:,1));
    st_y_1 = min(crds(:,2));
    st_y_2 = max(crds(:,2));
    len = st_y_2-st_y_1;
    no_div=round(len/w_tow);
    start_x = ones(no_div,1)*st_x;
    start_y = linspace(st_y_1,st_y_2,no_div);
elseif av_ang==90
    sgnV=sign(mean(mean(V)));
    sgnU=sgnV;
    % Start from bottom edge
    st_y = min(crds(:,2));
    st_x_1 = min(crds(:,1));
    st_x_2 = max(crds(:,1));
    len = st_x_2-st_x_1;
    no_div=round(len/w_tow);
    start_y = ones(no_div,1)*st_y;
    start_x = linspace(st_x_1,st_x_2,no_div);
end

towpath_1 = stream2(X,Y,U*sgnU,V*sgnV,start_x,start_y);
h = streamline(towpath_1);
set(h,'Color','blue');



% Plot the streamlines from the opposite edge
if av_ang==0
    sgnU=sign(mean(mean(U)));
    sgnV=sgnU;
    % Start from left edge
    st_x = max(crds(:,1));
    st_y_1 = min(crds(:,2));
    st_y_2 = max(crds(:,2));
    len = st_y_2-st_y_1;
    no_div=round(len/w_tow);
    start_x = ones(no_div,1)*st_x;
    start_y = linspace(st_y_1,st_y_2,no_div);
elseif av_ang==90
    sgnV=sign(mean(mean(V)));
    sgnU=sgnV;
    % Start from bottom edge
    st_y = max(crds(:,2));
    st_x_1 = min(crds(:,1));
    st_x_2 = max(crds(:,1));
    len = st_x_2-st_x_1;
    no_div=round(len/w_tow);
    start_y = ones(no_div,1)*st_y;
    start_x = linspace(st_x_1,st_x_2,no_div);
end

towpath_2 = stream2(X,Y,-U*sgnU,-V*sgnV,start_x,start_y);
i =streamline(towpath_2);
set(i,'Color','red');




title('stream lines')
axis equal







 % Now see the meaning of these potentials?

%   If you want the streamfunction only, use
%  psi = flowfun(u,v,'-')
%  (or psi = flowfun(u,v,'psi') or psi = flowfun(u,v,'stream') ).