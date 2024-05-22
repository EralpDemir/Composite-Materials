function [c,ceq,gradc,gradceq]= OPT_constraint_con(y,typLAM) 




switch typLAM
        
        case 'ORT'
            
            c = 2*y(1)^2-y(2)-1;
            

            
            dg_dy1 = 4*y(1);
            
            dg_dy2 = -1;
            
            
            gradc = [   dg_dy1;
                        dg_dy2 ];
            
            
        case 'GEN'    
 


            c=[ 2*y(1)^2*(1-y(3)) + 2*y(2)^2*(1+y(3)) + y(3)^2 + y(4)^2 - 4*y(1)*y(2)*y(4) - 1;
                y(1)^2 + y(2)^2 - 1;];






            dg1_dy1 = 4*y(1)*(1-y(3)) - 4*y(2)*y(4);
            dg2_dy1 = 2*y(1);

            dg1_dy2 = 4*y(2)*(1+y(3)) - 4*y(1)*y(4);
            dg2_dy2 = 2*y(2);

            dg1_dy3 = -2*y(1)^2 + 2*y(2)^2 +  2*y(3);
            dg2_dy3 = 0;

            dg1_dy4 = 2*y(4) - 4*y(1)*y(2);
            dg2_dy4 = 0;

            gradc =[    dg1_dy1, dg2_dy1;
                        dg1_dy2, dg2_dy2;
                        dg1_dy3, dg2_dy3;
                        dg1_dy4, dg2_dy4;];

                    
end



    



%         
ceq=[];

gradceq=[];


return
