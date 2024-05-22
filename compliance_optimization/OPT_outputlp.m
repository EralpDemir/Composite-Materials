function OPT_outputlp(nx, ny, L_x, ratio, el_vec_prop)

%OPT_ColPlot    plots the colorplot of elemental properties
%
% INPUTS:
% nx: number of elements in x-direction
% ny: number of elements in y-direction
% L_x: size of the laminate in x-direction
% ratio: length-to-width ratio
% el_vec_prop: elemental vector property

switch ratio
    
    case 1
        www = 13; % for ratio 1
        hhh = 12.4; % for ratio 1
        figure('Units','centimeters','Position',[0 0 www hhh])
        axis off
        axis equal % for ratio 1
        set(gcf,'paperunits','centimeters','papersize',[www hhh],...
            'paperposition',[0 0 www hhh])
        positionVector1 = [0.9, 1.85, 11.85, 10.2]; % for ratio 1
        
    case 3
        www = 26.2; % for ratio 3
        hhh = 9.7; % for ratio 3
        figure('Units','centimeters','Position',[0 0 www hhh])
        axis off
        set(gcf,'paperunits','centimeters','papersize',[www hhh],...
            'paperposition',[0 0 www hhh])
        positionVector1 = [0.9, 1.85, 25.1, 7.5]; % for ratio 3
end

axes('units','centimeters','position',positionVector1);
el_vec_prop_rshpd = reshape(el_vec_prop,[nx,ny])';
el_vec_prop_flpd = flip(el_vec_prop_rshpd);
x = linspace(0,L_x,nx);
y = linspace(0,L_x/ratio,ny);
imagesc(x,y,el_vec_prop_flpd)
colormap(jet(256))

%%%%%%%%%%%%%%%% Plotting mesh lines %%%%%%%%%%%%%%%%%%%%%%%%%%
dx = L_x/(nx-1); dy = L_x/ratio/(ny-1);
for i = -dx/2:dx:L_x+dx/2
    line([i i], [-dy/2 L_x/ratio+dy/2], 'Color','k');
end
for i = -dy/2:dy:L_x/ratio+dy/2
    line([-dx/2 L_x+dy/2], [i i], 'Color','k');
end
set(gca,'xtick',[],'ytick',[]);

end
