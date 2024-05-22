function [c,ceq,gradc,gradceq]= OPT_constraint_var(x,numel,typLAM) 




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
 
%         
ceq=[];

gradceq=[];


return
