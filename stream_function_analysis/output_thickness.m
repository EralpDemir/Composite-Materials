function output_thickness(crds,totnod,thick)

figure
hold on
    
for i = 1:totnod
    

    

    x0 = crds(i,1);
    y0 = crds(i,2);
    z0 = thick(i);

    
    plot3(x0, y0, z0, 'bo')
    

end

%axis equal

