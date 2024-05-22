function FEM_output(d,dof,conn,crds,numel,nnpe,eltyp,scale)

% Output the displacements
figure
hold on

% Division of edges
inc = 0.2;

switch eltyp
    case 'TRI'
        
        
        for iele=1:numel
            
            % Node numbering of the element
            nodes=conn(iele,:);
            
            clear del
            % Undeformed coordinates
            for j=1:1:nnpe
                x0(j)=crds(nodes(j),1);
                y0(j)=crds(nodes(j),2);
                del(dof*j-dof+1:1:dof*j,1)=d(dof*nodes(j)-dof+1:1:dof*nodes(j));
            end
            
            % deformed
            x1 = crds(conn(iele,1),1) + del(1)*scale;
            y1 = crds(conn(iele,1),2)+ del(2)*scale;
            x2 = crds(conn(iele,2),1) + del(3)*scale;
            y2 = crds(conn(iele,2),2) + del(4)*scale;
            x3 = crds(conn(iele,3),1) + del(5)*scale;
            y3 = crds(conn(iele,3),2) + del(6)*scale;
            
%             % undeformed
%             x1 = crds(conn(iele,1),1);
%             y1 = crds(conn(iele,1),2);
%             x2 = crds(conn(iele,2),1);
%             y2 = crds(conn(iele,2),2);
%             x3 = crds(conn(iele,3),1);
%             y3 = crds(conn(iele,3),2);
            
            
            for k=0:inc:1
                
                
                %edge 1-2
                ksi=k;
                eta=1-k;
                
                [N,~,area] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                
                L1=N(1,1)/area;L2=N(1,3)/area;L3=N(1,5)/area;
                
                
                x = L1*x1+L2*x2+L3*x3;
                
                y = L1*y1+L2*y2+L3*y3;
                
                % plot3(x,y,z,'r.')
                plot(x , y, 'k.')
                
                %edge 2-3
                ksi=0;
                eta=k;
                
                [N,~,area] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                
                L1=N(1,1)/area;L2=N(1,3)/area;L3=N(1,5)/area;
                
                
                x = L1*x1+L2*x2+L3*x3;
                
                y = L1*y1+L2*y2+L3*y3;
                
                % plot3(x,y,z,'r.')
                plot(x , y, 'k.')
                
                %edge 3-1
                ksi=1-k;
                eta=0;
                
                [N,~,area] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                
                L1=N(1,1)/area;L2=N(1,3)/area;L3=N(1,5)/area;
                
                
                x = L1*x1+L2*x2+L3*x3;
                
                y = L1*y1+L2*y2+L3*y3;
                
                % plot3(x,y,z,'r.')
                plot(x , y, 'k.')
                
            end
            
        end
        
        
    case 'HCT'
        
        %     function FEM_output_hct(d,numel,conn,crds,dofst,scale)
        
        % % Output the displacements
        %figure
        % hold on
        
        
        % For each element
        for iele=1:numel
            
            % compute the element bending stiffness matrix and rhs
            j1=conn(iele,1); j2=conn(iele,2); j3=conn(iele,3);   % global element labels
            j4=conn(iele,4); j5=conn(iele,5); j6=conn(iele,6);
            
            % position in the global matrix
            e(1) = dofst(j1);
            e(2) = dofst(j1)+1;
            e(3) = dofst(j1)+2;
            e(4) = dofst(j2);
            e(5) = dofst(j2)+1;
            e(6) = dofst(j2)+2;
            e(7) = dofst(j3);
            e(8) = dofst(j3)+1;
            e(9) = dofst(j3)+2;
            e(10)= dofst(j4);
            e(11)= dofst(j5);
            e(12)= dofst(j6);
            
            del=d(e);
            
            % Element coordinates
            x1=crds(j1,1);
            y1=crds(j1,2);
            x2=crds(j2,1);
            y2=crds(j2,2);
            x3=crds(j3,1);
            y3=crds(j3,2);
            
            %----
            % interior node
            %-----
            x7=(x1+x2+x3)/3.0; y7=(y1+y2+y3)/3.0;
            
            
            [c] = FEM_hctsys (x1,y1, x2,y2, x3,y3, del);
            
            
            % For each sub-element
            for isub=1:1:3
                
                for ksi=0:inc:1
                    eta=0;
                %for ksi=0:0.5:1
                    
                    %for eta=0:0.5:(1-ksi)
                        
                        
                        if(isub==1)       % triangle 1-2-7
                            
                            x = x1 + (x2-x1)*ksi + (x7-x1)*eta;
                            y = y1 + (y2-y1)*ksi + (y7-y1)*eta;
                            
                            
                            
                            psi = c(1) + c(2)*x   +c(3)*y ...
                                + c(4)*x^2 +c(5)*x*y   + c(6)*y^2 ...
                                + c(7)*x^3 +c(8)*x^2*y+c(9)*x*y^2 +c(10)*y^3;
                            
                            
                            
                        elseif(isub==2)  % triangle 2-3-7
                            
                            x = x2 + (x3-x2)*ksi + (x7-x2)*eta;
                            y = y2 + (y3-y2)*ksi + (y7-y2)*eta;
                            
                            
                            
                            psi = c(11) + c(12)*x   +c(13)*y ...
                                + c(14)*x^2 +c(15)*x*y   + c(16)*y^2 ...
                                + c(17)*x^3 +c(18)*x^2*y+c(19)*x*y^2 +c(20)*y^3;
                            
                            
                            
                        else            % triangle 3-1-7
                            
                            x = x3 + (x1-x3)*ksi + (x7-x3)*eta;
                            y = y3 + (y1-y3)*ksi + (y7-y3)*eta;
                            
                            
                            
                            psi = c(21) + c(22)*x   +c(23)*y ...
                                + c(24)*x^2 +c(25)*x*y   + c(26)*y^2 ...
                                + c(27)*x^3 +c(28)*x^2*y+c(29)*x*y^2 +c(30)*y^3;
                            
                            
                        end
                        
                        w = scale * psi;
                        
                        
                        plot3(x,y,w,'k.');
                        
                    %end
                    
                    
                %end
                end
            
            end
            
        end
        
        
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
        
        
        
        
    case 'BFS+LAG'
        
        
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
                uvw=scale*N*del;
                
                %        plot3(x,y,z,'r.')
                plot3(x+uvw(1),y+uvw(2),uvw(3),'k.')
                
                %
                % Edge 2-3
                eta=k;ksi=1;
                x=a*ksi+xc;y=b*eta+yc;
                
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uvw=scale*N*del;
                %        plot3(x,y,z,'r.')
               plot3(x+uvw(1),y+uvw(2),uvw(3),'k.')
                
                %
                % Edge 3-4
                eta=1;ksi=k;
                x=a*ksi+xc;y=b*eta+yc;
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uvw=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot3(x+uvw(1),y+uvw(2),uvw(3),'k.')
                
                %
                % Edge 4-1
                eta=k;ksi=-1;
                x=a*ksi+xc;y=b*eta+yc;
                
                [N] = FEM_nmat(ksi,eta,crds,conn,eltyp,iele);
                uvw=scale*N*del;
                %        plot3(x,y,z,'r.')
                plot3(x+uvw(1),y+uvw(2),uvw(3),'k.')
                
            end
            
        end        
        
end

axis equal

end





