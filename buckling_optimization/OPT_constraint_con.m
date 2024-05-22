function [c,ceq,gradc,gradceq]= OPT_constraint_con(x) 



% option-1
t=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.50, 0.75, 1];

% % option-2
% t=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1];

nt=size(t,2);
ncon=2*nt+1;



c=[];









 %============================ First Formula============================ 
 g0 = 5*(x(1) - x(3))^2 - 2*(1+x(2)-2*x(1)^2);



 dg_dy1 = 18*x(1) - 10*x(3);

 dg_dy2 = -2;

 dg_dy3 = -10*x(1) + 10*x(3);

 dg_dy4 = 0;
 %============================

c = g0;


 gradc = [  dg_dy1;
            dg_dy2;
            dg_dy3;
            dg_dy4];



 %============================ 2nd and 3rd Formula============================
 for i=1:1:nt



      % For the first constraint
     g1 =(x(2)-4*t(i)*x(1)+1+2*t(i)^2)^3 - 4*(1+2*abs(t(i))+t(i)^2)^2*(x(4)-4*t(i)*x(3)+1+2*t(i)^2);

     dg_dy1 = -12*t(i)*(x(2)-4*t(i)*x(1)+1+2*t(i)^2)^2;

     dg_dy2 = 3*(x(2)-4*t(i)*x(1)+1+2*t(i)^2)^2;

     dg_dy3 = 16*t(i)*(1+2*abs(t(i))+t(i)^2)^2;

     dg_dy4 = -4*(1+2*abs(t(i))+t(i)^2)^2;



     c= [c; g1;];


     dg = [ dg_dy1;
            dg_dy2;
            dg_dy3;
            dg_dy4];

     gradc = [gradc, dg];


     % For the second constraint
     g2 =(4*t(i)*x(1)-x(2)+1+4*abs(t(i)))^3 - 4*(1+2*abs(t(i))+t(i)^2)^2*(4*t(i)*x(3) - x(4) + 1 + 4*abs(t(i)));

     dg_dy1 = 12*t(i)*(4*t(i)*x(1)-x(2)+1+4*abs(t(i)))^2 ;

     dg_dy2 = -3*(4*t(i)*x(1)-x(2)+1+4*abs(t(i)))^2;

     dg_dy3 = -16*t(i)*(1+2*abs(t(i))+t(i)^2)^2;

     dg_dy4 = 4*(1+2*abs(t(i))+t(i)^2)^2;




     c= [c; g2;];

     dg = [ dg_dy1;
            dg_dy2;
            dg_dy3;
            dg_dy4];

     gradc = [gradc, dg];





 end 



  

    
ceq = [];

gradceq = [];

return
