function [dndx,dndy,detj] = jacobian(nnpe,nqpt,dN_qpt,elcrds)
% LSC_sfder computes shape function derivatives 

detj = zeros(1,nqpt);
dndx = zeros(nqpt,nnpe);
dndy = zeros(nqpt,nnpe);

dn_dksi(1:nqpt,1:nnpe) = dN_qpt(1:nqpt,1,1:nnpe);
dn_deta(1:nqpt,1:nnpe) = dN_qpt(1:nqpt,2,1:nnpe);

x = elcrds(:,1);
y = elcrds(:,2);

for i = 1:1:nqpt
    
    tjac = zeros(2,2);
         
    for j = 1:1:nnpe  
	  
        tjac(1,1) = tjac(1,1) + dn_dksi(i,j)*x(j);                           
        tjac(1,2) = tjac(1,2) + dn_deta(i,j)*x(j);	  
                         
        tjac(2,1) = tjac(2,1) + dn_dksi(i,j)*y(j);                          
        tjac(2,2) = tjac(2,2) + dn_deta(i,j)*y(j);                         
 	  
    end                                                         
                                                                    
    detj(i) = det(tjac);

    tjaci = inv(tjac);                                      
                   
    for  j = 1:1:nnpe
	
        dndx(i,j) = tjaci(1,1)*dn_dksi(i,j)+tjaci(2,1)*dn_deta(i,j);            
        dndy(i,j) = tjaci(1,2)*dn_dksi(i,j)+tjaci(2,2)*dn_deta(i,j);       

    end      
  
end
 

