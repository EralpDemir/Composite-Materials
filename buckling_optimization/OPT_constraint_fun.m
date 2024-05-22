function [c]= OPT_constraint_fun(x) 



% option-1
t=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.50, 0.75, 1];

% % option-2
% t=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1];

nt=size(t,2);










 %============================ First Formula============================ 
 g0 = 5*(x(1) - x(3))^2 - 2*(1+x(2)-2*x(1)^2);


c = g0;



 %============================ 2nd and 3rd Formula============================
 for i=1:1:nt



      % For the first constraint
     g1 =(x(2)-4*t(i)*x(1)+1+2*t(i)^2)^3 - 4*(1+2*abs(t(i))+t(i)^2)^2*(x(4)-4*t(i)*x(3)+1+2*t(i)^2);



     c= [c; g1;];





     % For the second constraint
     g2 =(4*t(i)*x(1)-x(2)+1+4*abs(t(i)))^3 - 4*(1+2*abs(t(i))+t(i)^2)^2*(4*t(i)*x(3) - x(4) + 1 + 4*abs(t(i)));






     c= [c; g2;];






 end 



  



return
