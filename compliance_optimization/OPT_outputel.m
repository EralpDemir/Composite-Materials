function OPT_outputel(conn,crds,numel,nnpe,angPly,izPly)

figure
hold on
    
for i = 1:numel
    
    % Node numbering of the element
    nodes = conn(i,:);
    
    
    % Undeformed coordinates
    for j = 1:1:nnpe
        x0(j) = crds(nodes(j),1);
        y0(j) = crds(nodes(j),2);
        %         del(6*j,1)=d(6*nodes(j));
        %         del(6*j-1,1)=d(6*nodes(j)-1);
        %         del(6*j-2,1)=d(6*nodes(j)-2);
        %         del(6*j-3,1)=d(6*nodes(j)-3);
        %         del(6*j-4,1)=d(6*nodes(j)-4);
        %         del(6*j-5,1)=d(6*nodes(j)-5);
    end
    
    % Plot the undeformed geometry of elements
    line([x0(1:nnpe) x0(1)],[y0(1:nnpe) y0(1)],'Color','k');
    
    xc = mean(x0);
    yc = mean(y0);
    
    dxy = 0.25*sqrt((max(y0)-min(y0))^2+(max(x0)-min(x0))^2);
    
    ang = angPly(i,izPly)*pi/180;
    
    % Plot the fiber angles
    
    % Find the start & end coordinates of fibers
    xs = xc+dxy*cos(pi+ang);
    ys = yc+dxy*sin(pi+ang);
    xe = xc+dxy*cos(ang);
    ye = yc+dxy*sin(ang);
    
    plot([xs xe], [ys ye], 'r')
    

end

axis equal

