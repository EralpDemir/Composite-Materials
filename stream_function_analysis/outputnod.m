function outputnod(crds,conn,nnpe,numel,angPlyNodes)


% angles are given in radians


% Output the displacements
% Plot after output
figure
hold on



dxy_=zeros(1,numel);

for i=1:1:numel

  % Node numbering of the element
    nodes=conn(i,:);
    
    
    
        % Undeformed coordinates
    for j=1:1:nnpe
        x0(j)=crds(nodes(j),1);
        y0(j)=crds(nodes(j),2);
%         del(6*j,1)=d(6*nodes(j));
%         del(6*j-1,1)=d(6*nodes(j)-1);
%         del(6*j-2,1)=d(6*nodes(j)-2);
%         del(6*j-3,1)=d(6*nodes(j)-3);
%         del(6*j-4,1)=d(6*nodes(j)-4);
%         del(6*j-5,1)=d(6*nodes(j)-5);
    end
    
    % Plot the undeformed geometry of elements
    line([x0(1:nnpe) x0(1)],[y0(1:nnpe) y0(1)],'Color','k');


    dxy_(i)=0.25*sqrt((max(y0)-min(y0))^2+(max(x0)-min(x0))^2);
    
end
    
    
dxy=mean(dxy_);


for i=1:1:size(crds,1)
    
   
    xc=crds(i,1);
    yc=crds(i,2);
    
   
  
    
    
    

    % Read the inputs

    angply=angPlyNodes(i);
    
    
    % Plot the fiber angles 
    
    % Find the start & end coordinates of fibers
    xs=xc+dxy*cos(angply+pi);
    ys=yc+dxy*sin(angply+pi);
    xe=xc+dxy*cos(angply);
    ye=yc+dxy*sin(angply);

    
    plot([xs xe], [ys ye], 'b')
    
    
    
end

axis equal

