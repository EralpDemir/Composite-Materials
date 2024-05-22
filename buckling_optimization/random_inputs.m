% Random inputs

function [All_xi_A0, All_xi_D0] = random_inputs(numel)


All_xi_A0=zeros(4,numel);

All_xi_D0=zeros(4,numel);


for i=1:1:numel
    
   flag = true;
   
   while flag
       
       xi_A1=(2*rand-1);
       
       xi_A3=(2*rand-1);
       
       xi_D1=(2*rand-1);
       
       xi_D3=(2*rand-1);
    

       x = [xi_A1, xi_A3, xi_D1, xi_D3];
       
       [c,~] = OPT_constraint_con(x);
       
       flag=false;
       for j=1:1:19
        if c(j)>0
            flag=true;
        end
       end
       
           
       
   end 
   
   All_xi_A0(1,i) = xi_A1;
   
   All_xi_A0(3,i) = xi_A3;
   
   All_xi_D0(1,i) = xi_D1;
   
   All_xi_D0(3,i) = xi_D3;
   
   
    
end