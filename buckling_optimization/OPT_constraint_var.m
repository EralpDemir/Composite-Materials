function [c,ceq,gradc,gradceq]= OPT_constraint_var(x,numel) 



% option-1
t=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.50, 0.75, 1];

% 
% % option-2
% t=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1];


nt=size(t,2);
ncon=2*nt+1;

c=zeros(ncon*numel,1);



gradc=zeros(4*numel,ncon*numel);



for iele=1:1:numel

    y=zeros(4,1);

     y(1:4,1)=x(4*iele-3:1:4*iele,1);
     
     
     g = [];


     %============================ First Formula============================ 
     g0 = 5*(y(1) - y(3))^2 - 2*(1+y(2)-2*y(1)^2);



     dg_dy1 = 18*y(1) - 10*y(3);

     dg_dy2 = -2;

     dg_dy3 = -10*y(1) + 10*y(3);

     dg_dy4 = 0;
     %============================

    %c((iele-1)*ncon+1, 1) = g0;

    g = g0;
    

     gradg = [  dg_dy1;
                dg_dy2;
                dg_dy3;
                dg_dy4];



     %============================ 2nd and 3rd Formula============================
     for i=1:1:nt

      

         % For the first constraint
         g1 =(y(2)-4*t(i)*y(1)+1+2*t(i)^2)^3 - 4*(1+2*abs(t(i))+t(i)^2)^2*(y(4)-4*t(i)*y(3)+1+2*t(i)^2);

         dg_dy1 = -12*t(i)*(y(2)-4*t(i)*y(1)+1+2*t(i)^2)^2;

         dg_dy2 = 3*(y(2)-4*t(i)*y(1)+1+2*t(i)^2)^2;

         dg_dy3 = 16*t(i)*(1+2*abs(t(i))+t(i)^2)^2;

         dg_dy4 = -4*(1+2*abs(t(i))+t(i)^2)^2;


%          ind1 = 2*i-1;

%          c((iele-1)*ncon+ind1+1,1) = g1;

         g = [  g;
                g1;];


         dg = [ dg_dy1;
                dg_dy2;
                dg_dy3;
                dg_dy4];

         gradg = [gradg, dg];


         % For the second constraint
         g2 =(4*t(i)*y(1)-y(2)+1+4*abs(t(i)))^3 - 4*(1+2*abs(t(i))+t(i)^2)^2*(4*t(i)*y(3) - y(4) + 1 + 4*abs(t(i)));

         dg_dy1 = 12*t(i)*(4*t(i)*y(1)-y(2)+1+4*abs(t(i)))^2 ;

         dg_dy2 = -3*(4*t(i)*y(1)-y(2)+1+4*abs(t(i)))^2;

         dg_dy3 = -16*t(i)*(1+2*abs(t(i))+t(i)^2)^2;

         dg_dy4 = 4*(1+2*abs(t(i))+t(i)^2)^2;


%          ind2 = 2*i;

%          c((iele-1)*ncon+ind2+1,1) = g2;

         g = [g;
              g2;];  

         dg = [ dg_dy1;
                dg_dy2;
                dg_dy3;
                dg_dy4];

         gradg = [gradg, dg];





     end 

     
     
    c(ncon*(iele-1)+1:1:ncon*iele,1)=g; 
     
     

    gradc(4*(iele-1)+1:1:4*iele,(iele-1)*ncon+1:1:iele*ncon)=gradg;

end

%gradc =zeros(4*numel,ncon*numel);


%         
ceq=[];

gradceq=[];


return