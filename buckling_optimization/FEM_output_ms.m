function FEM_output_ms(d,dof,conn,crds,numel,nnpe,eltyp,scale)

% Output the displacements
figure
hold on

% Division of edges
inc = 1;

switch eltyp
 
      
        
        
    case 'REC'
        
        
        for iele=1:numel
            
            a = abs(crds(conn(iele,2),1)-crds(conn(iele,1),1))/2;
            b = abs(crds(conn(iele,4),2)-crds(conn(iele,1),2))/2;
            
            xc = mean([  crds(conn(iele,1),1) crds(conn(iele,2),1) ...
                crds(conn(iele,3),1) crds(conn(iele,4),1)]);
            yc = mean([  crds(conn(iele,1),2) crds(conn(iele,2),2) ...
                crds(conn(iele,3),2) crds(conn(iele,4),2)]);
            
            % Node numbering of the element
            nodes=conn(iele,:);
            
            clear del
            % Undeformed coordinates
            for j=1:1:nnpe
                x0(j)=crds(nodes(j),1);
                y0(j)=crds(nodes(j),2);
                del(dof*j-dof+1:1:dof*j,1)=d(dof*nodes(j)-dof+1:1:dof*nodes(j));
            end
            
            line([x0(1:nnpe) x0(1)],[y0(1:nnpe) y0(1)],'Color','k');
            
            % for all ksi and eta at the edges
            for k=-1:inc:1
                
                % Edge 1-2
                eta=-1;ksi=k;
                x=a*ksi+xc;y=b*eta+yc;
                
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uv=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot(x+uv(1),y+uv(2),'k.')
                
                %
                % Edge 2-3
                eta=k;ksi=1;
                x=a*ksi+xc;y=b*eta+yc;
                
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uv=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot(x+uv(1),y+uv(2),'k.')
                
                %
                % Edge 3-4
                eta=1;ksi=k;
                x=a*ksi+xc;y=b*eta+yc;
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uv=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot(x+uv(1),y+uv(2),'k.')
                
                %
                % Edge 4-1
                eta=k;ksi=-1;
                x=a*ksi+xc;y=b*eta+yc;
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uv=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot(x+uv(1),y+uv(2),'k.')
                
            end
            
        end
        
    case 'BFS'
        
        
        for iele=1:numel
            
            a = abs(crds(conn(iele,2),1)-crds(conn(iele,1),1))/2;
            b = abs(crds(conn(iele,4),2)-crds(conn(iele,1),2))/2;
            
            xc = mean([  crds(conn(iele,1),1) crds(conn(iele,2),1) ...
                crds(conn(iele,3),1) crds(conn(iele,4),1)]);
            yc = mean([  crds(conn(iele,1),2) crds(conn(iele,2),2) ...
                crds(conn(iele,3),2) crds(conn(iele,4),2)]);
            
            % Node numbering of the element
            nodes=conn(iele,:);
            
            clear del
            % Undeformed coordinates
            for j=1:1:nnpe
                x0(j)=crds(nodes(j),1);
                y0(j)=crds(nodes(j),2);
                del(dof*j-dof+1:1:dof*j,1)=d(dof*nodes(j)-dof+1:1:dof*nodes(j));
            end
            
            line([x0(1:nnpe) x0(1)],[y0(1:nnpe) y0(1)],'Color','k');
            
            % for all ksi and eta at the edges
            for k=-1:inc:1
                
                % Edge 1-2
                eta=-1;ksi=k;
                x=a*ksi+xc;y=b*eta+yc;
                
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                w=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot3(x,y,w,'k.')
                
                %
                % Edge 2-3
                eta=k;ksi=1;
                x=a*ksi+xc;y=b*eta+yc;
                
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                w=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot3(x,y,w,'k.')
                
                %
                % Edge 3-4
                eta=1;ksi=k;
                x=a*ksi+xc;y=b*eta+yc;
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                w=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot3(x,y,w,'k.')
                
                %
                % Edge 4-1
                eta=k;ksi=-1;
                x=a*ksi+xc;y=b*eta+yc;
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                w=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot3(x,y,w,'k.')
                
            end
            
        end
        
        
        
        
end

axis equal

end





