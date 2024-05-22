function [c,ceq,gradc,gradceq]= OPT_constraint_var(x,typOPT,numel,typLAM) 




        


if typOPT==1
    
     switch typLAM
        
        case 'ORT'
            
            c=zeros(numel,1);

            gradc=zeros(2*numel,numel);
                        
            
            
            for iele=1:1:numel
            
                 y(1:2,1)=x(2*iele-1:1:2*iele,1);
                
                 c(iele) = 2*y(1)^2-y(2)-1;
            

            
                 dg_dy1 = 4*y(1);
            
                 dg_dy2 = -1;
                 
                 gradg = [  dg_dy1;
                            dg_dy2 ];
                
                        
                gradc(2*iele-1:1:2*iele,iele)=gradg;
            
            end
            
            
        case 'GEN'
    
            c=zeros(2*numel,1);

            gradc=zeros(4*numel,2*numel);
            
            
            for iele=1:1:numel

                y(1:4,1)=x(4*iele-3:1:4*iele,1);

                g=[ 2*y(1)^2*(1-y(3)) + 2*y(2)^2*(1+y(3)) + y(3)^2 + y(4)^2 - 4*y(1)*y(2)*y(4) - 1;
                    y(1)^2 + y(2)^2 - 1;];

                c(2*iele-1:1:2*iele,1)=g;




                dg1_dy1 = 4*y(1)*(1-y(3)) - 4*y(2)*y(4);
                dg2_dy1 = 2*y(1);

                dg1_dy2 = 4*y(2)*(1+y(3)) - 4*y(1)*y(4);
                dg2_dy2 = 2*y(2);

                dg1_dy3 = -2*y(1)^2 + 2*y(2)^2 +  2*y(3);
                dg2_dy3 = 0;

                dg1_dy4 = 2*y(4) - 4*y(1)*y(2);
                dg2_dy4 = 0;

                gradg =[    dg1_dy1, dg2_dy1;
                            dg1_dy2, dg2_dy2;
                            dg1_dy3, dg2_dy3;
                            dg1_dy4, dg2_dy4;];


               gradc(4*iele-3:1:4*iele,2*iele-1:1:2*iele)=gradg;         

            end
     end
    
    
     
elseif typOPT==2

        
     switch typLAM
        
        case 'ORT'
            
            c=zeros(numel,1);

            gradc=zeros(2*numel,numel);
                        
            
            
            for iele=1:1:numel
            
                 y(1:2,1)=x(2*iele-1:1:2*iele,1);
                
                 g = 2*y(1)^2-y(2)-1;
            

                 c(iele)=g;
            
                 dg_dy1 = 4*y(1);
            
                 dg_dy2 = -1;
                 
                 gradg = [  dg_dy1;
                            dg_dy2 ];
                
                        
                gradc(2*iele-1:1:2*iele,iele)=gradg;
            
            end
            
            
        case 'GEN'
    
            c=zeros(2*numel,1);

            gradc=zeros(4*numel,2*numel);
            
    
            for iele=1:1:numel
                y(1:4,1)=x(4*iele-3:1:4*iele,1);


                g=[ 2*y(1)^2*(1-y(3)) + 2*y(2)^2*(1+y(3)) + y(3)^2 + y(4)^2 - 4*y(1)*y(2)*y(4) - 1;
                y(1)^2 + y(2)^2 - 1;];


                c(2*iele-1:1:2*iele,1)=g;

                dg1_dy1 = 4*y(1)*(1-y(3)) - 4*y(2)*y(4);
                dg2_dy1 = 2*y(1);

                dg1_dy2 = 4*y(2)*(1+y(3)) - 4*y(1)*y(4);
                dg2_dy2 = 2*y(2);

                dg1_dy3 = -2*y(1)^2 + 2*y(2)^2 +  2*y(3);
                dg2_dy3 = 0;

                dg1_dy4 = 2*y(4) - 4*y(1)*y(2);
                dg2_dy4 = 0;

                gradg =[    dg1_dy1, dg2_dy1;
                            dg1_dy2, dg2_dy2;
                            dg1_dy3, dg2_dy3;
                            dg1_dy4, dg2_dy4;];

                gradc(4*iele-3:1:4*iele,2*iele-1:1:2*iele)=gradg;

            end
     end
% elseif typOPT==3
%     
%     c=[ 2*y(1)^2*(1-y(3)) + 2*y(2)^2*(1+y(3)) + y(3)^2 + y(4)^2 - 4*y(1)*y(2)*y(4) - 1;
%         y(1)^2 + y(2)^2 - 1;
%         2*y(5)^2*(1-y(7)) + 2*y(6)^2*(1+y(7)) + y(7)^2 + y(8)^2 - 4*y(5)*y(6)*y(8 )- 1;
%         y(5)^2 + y(6)^2 - 1;];
%         
%     dg1_dy1 = 4*y(1)*(1-y(3)) - 4*y(2)*y(4);
%     dg2_dy1 = 2*y(1);
%     dg3_dy1 = 0;
%     dg4_dy1 = 0;
%     
%     dg1_dy2 = 4*y(2)*(1+y(3)) + 2*y(2)^2 - 4*y(1)*y(4);
%     dg2_dy2 = 2*y(2);
%     dg3_dy2 = 0;
%     dg4_dy2 = 0;
%     
%     dg1_dy3 = -2*y(1)^2 + 2*y(2)^2 +  2*y(3);
%     dg2_dy3 = 0;
%     dg3_dy3 = 0;
%     dg4_dy3 = 0;
%     
%     dg1_dy4 = 2*y(4) - 4*y(1)*y(2);
%     dg2_dy4 = 0;
%     dg3_dy4 = 0;
%     dg4_dy4 = 0;
%     
%     
%     
% 
%     dg1_dy5 = 0;
%     dg2_dy5 = 0;
%     dg3_dy5 = 4*y(5)*(1-y(7)) - 4*y(6)*y(8);
%     dg4_dy5 = 2*y(5);
% 
% 
%     dg1_dy6 = 0;
%     dg2_dy6 = 0;
%     dg3_dy6 = 4*y(6)*(1+y(7)) + 2*y(6)^2 - 4*y(5)*y(8);
%     dg4_dy6 = 2*y(6);
% 
% 
%     dg1_dy7 = 0;
%     dg2_dy7 = 0;
%     dg3_dy7 = -2*y(5)^2 + 2*y(6)^2 + 2*y(7);
%     dg4_dy7 = 0;
% 
% 
%     dg1_dy8 = 0;
%     dg2_dy8 = 0;
%     dg3_dy8 = 2*y(8) - 4*y(5)*y(6);
%     dg4_dy8 = 0;
%     
%     
%     
% 
%     gradc =[    dg1_dy1, dg2_dy1, dg3_dy1, dg4_dy1;
%                 dg1_dy2, dg2_dy2, dg3_dy2, dg4_dy2;
%                 dg1_dy3, dg2_dy3, dg3_dy3, dg4_dy3;
%                 dg1_dy4, dg2_dy4, dg3_dy4, dg4_dy4;
% 
% 
%                 dg1_dy5, dg2_dy5, dg3_dy5, dg4_dy5;
%                 dg1_dy6, dg2_dy6, dg3_dy6, dg4_dy6;
%                 dg1_dy7, dg2_dy7, dg3_dy7, dg4_dy7;
%                 dg1_dy8, dg2_dy8, dg3_dy8, dg4_dy8;
% 
%                 ];

end


% c=[ 2*y(1)^2*(1-y(3)) + 2*y(2)^2*(1+y(3)) + y(3)^2 + y(4)^2 - 4*y(1)*y(2)*y(4) - 1;
%     y(1)^2 + y(2)^2 - 1;
%     y(3) - 1;
%    -y(3) - 1;];



% dc_dy1 = 4*y(1);
% dc_dy2 = -1;

% dg2_dy1 = 1;
% dg3_dy1 = -1;
% dg4_dy1 = 0;
% dg5_dy1 = 0;

% dg2_dy2 = 0;
% dg3_dy2 = 0;
% dg4_dy2 = 1;
% dg5_dy2 = -1;






% dg1_dy1 = 4*y(1)*(1-y(3)) - 4*y(2)*y(4);
% dg2_dy1 = 2*y(1);
% dg3_dy1 = 0;
% dg4_dy1 = 0;

% dg1_dy2 = 4*y(2)*(1+y(3)) + 2*y(2)^2 - 4*y(1)*y(4);
% dg2_dy2 = 2*y(2);
% dg3_dy2 = 0;
% dg4_dy2 = 0;

% dg1_dy3 = -2*y(1)^2 + 2*y(2)^2 +  2*y(3);
% dg2_dy3 = 0;
% dg3_dy3 = 1;
% dg4_dy3 = -1;

% dg1_dy4 = 2*y(4) - 4*y(1)*y(2);
% dg2_dy4 = 0;
% dg3_dy4 = 0;
% dg4_dy4 = 0;


% gradc =[    dg1_dy1, dg2_dy1, dg3_dy1, dg4_dy1;
%             dg1_dy2, dg2_dy2, dg3_dy2, dg4_dy2;
%             dg1_dy3, dg2_dy3, dg3_dy3, dg4_dy3;
%             dg1_dy4, dg2_dy4, dg3_dy4, dg4_dy4;
%             
%             zeros(8,4);
%             
%             dg1_dy5, dg2_dy5, dg3_dy5, dg4_dy5;
%             dg1_dy6, dg2_dy6, dg3_dy6, dg4_dy6;
%             dg1_dy7, dg2_dy7, dg3_dy7, dg4_dy7;
%             dg1_dy8, dg2_dy8, dg3_dy8, dg4_dy8;
%             
%             ];
%         
% 
% gradc =[    dg1_dy1, dg2_dy1, dg3_dy1, dg4_dy1;
%             dg1_dy2, dg2_dy2, dg3_dy2, dg4_dy2;
%             dg1_dy3, dg2_dy3, dg3_dy3, dg4_dy3;
%             dg1_dy4, dg2_dy4, dg3_dy4, dg4_dy4;
%          
%             
%             ];        
%     
%          
%             
%             ];  

% 
% gradc=[ dc_dy1; 
%         dc_dy2];


% gradc =[    dg1_dy1, dg2_dy1, dg3_dy1, dg4_dy1, dg5_dy1;
%             dg1_dy2, dg2_dy2, dg3_dy2, dg4_dy2, dg5_dy2;];
%         
ceq=[];

gradceq=[];