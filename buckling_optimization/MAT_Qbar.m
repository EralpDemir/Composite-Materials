
function [Q_bar] = MAT_Qbar(nPly, angPly, Q11, Q12, Q22, Q66)

conv=pi/180;

Q_bar=zeros(3,3,nPly);


%Q_bar=zeros(3,3,nPly);

for i=1:1:nPly
    
    % Conversion to radians
    theta=angPly(i)*conv;
    
%CALCULATE REDUCED STIFFNESS MATRIX IN ORDER TO COMPUTE STRESS FOR EACH
%LAYER   ZIGMA=Q_BAR*(ex0+zkapax)

Q_bar(1,1,i) = (cos(theta)^4)*Q11 + (sin(theta)^4)*Q22 + (cos(theta)^2)*(sin(theta)^2)*(2*Q12+4*Q66);
Q_bar(1,2,i) = (cos(theta)^2)*(sin(theta)^2)*(Q11 + Q22 - 4*Q66) + (cos(theta)^4 + sin(theta)^4)* Q12;
Q_bar(1,3,i) = (cos(theta)^3)*sin(theta)*(Q11-Q12-2*Q66) + cos(theta)*(sin(theta)^3)*(Q12-Q22+2*Q66);
Q_bar(2,2,i) = (sin(theta)^4)*Q11+ (cos(theta)^2)*(sin(theta)^2)*(2*Q12 + 4*Q66) + (cos(theta)^4)*Q22;
Q_bar(2,3,i) = (sin(theta)^3)*cos(theta)*(Q11-Q12-2*Q66) + (cos(theta)^3)*sin(theta)*(Q12-Q22+2*Q66);
Q_bar(3,3,i) = (cos(theta)^2)*(sin(theta)^2)*(Q11 + Q22-2*Q12-2*Q66) + (sin(theta)^4 + cos(theta)^4)*Q66;

Q_bar(2,1,i)=Q_bar(1,2,i);
Q_bar(3,2,i)=Q_bar(2,3,i);
Q_bar(3,1,i)=Q_bar(1,3,i);

end
