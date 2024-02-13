function [r,y] = beddington (x,rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix,Imatrix,Rmatrix,Fmatrix)

r(1) =  rcons(:,1) + delta*rrain(:,1)*rain2 + delta*rtemp(:,1)*temp2 - Imatrix*alpha1(:,1) - Rmatrix*alpha2(:,1) - Fmatrix*alpha3(:,1) - Pmatrix.*Imatrix*tau1(:,1) - Pmatrix.*Rmatrix*tau2(:,1) - Pmatrix.*Fmatrix*tau3(:,1);
r(2) =  rcons(:,2) + delta*rrain(:,2)*rain2 + delta*rtemp(:,2)*temp2 - Imatrix*alpha1(:,2) - Rmatrix*alpha2(:,2) - Fmatrix*alpha3(:,2) - Pmatrix.*Imatrix*tau1(:,2) - Pmatrix.*Rmatrix*tau2(:,2) - Pmatrix.*Fmatrix*tau3(:,2);
r(3) =  rcons(:,3) + delta*rrain(:,3)*rain2 + delta*rtemp(:,3)*temp2 - Imatrix*alpha1(:,3) - Rmatrix*alpha2(:,3) - Fmatrix*alpha3(:,3) - Pmatrix.*Imatrix*tau1(:,3) - Pmatrix.*Rmatrix*tau2(:,3) - Pmatrix.*Fmatrix*tau3(:,3);


y(1) = x(1) + rcons(:,1) + delta*rrain(:,1)*rain2 + delta*rtemp(:,1)*temp2 - Imatrix*alpha1(:,1) - Rmatrix*alpha2(:,1) - Fmatrix*alpha3(:,1) - Pmatrix.*Imatrix*tau1(:,1) - Pmatrix.*Rmatrix*tau2(:,1) - Pmatrix.*Fmatrix*tau3(:,1);
y(2) = x(2) + rcons(:,2) + delta*rrain(:,2)*rain2 + delta*rtemp(:,2)*temp2 - Imatrix*alpha1(:,2) - Rmatrix*alpha2(:,2) - Fmatrix*alpha3(:,2) - Pmatrix.*Imatrix*tau1(:,2) - Pmatrix.*Rmatrix*tau2(:,2) - Pmatrix.*Fmatrix*tau3(:,2);
y(3) = x(3) + rcons(:,3) + delta*rrain(:,3)*rain2 + delta*rtemp(:,3)*temp2 - Imatrix*alpha1(:,3) - Rmatrix*alpha2(:,3) - Fmatrix*alpha3(:,3) - Pmatrix.*Imatrix*tau1(:,3) - Pmatrix.*Rmatrix*tau2(:,3) - Pmatrix.*Fmatrix*tau3(:,3);
% 
% r(1) =  rcons(:,1) + delta*rrain(:,1)*rain2 + delta*rtemp(:,1)*temp2 - Imatrix*alpha1(:,1) - Fmatrix*alpha2(:,1) - Rmatrix*alpha3(:,1) - Pmatrix.*Imatrix*tau1(:,1) - Pmatrix.*Fmatrix*tau2(:,1) - Pmatrix.*Rmatrix*tau3(:,1);
% r(2) =  rcons(:,2) + delta*rrain(:,2)*rain2 + delta*rtemp(:,2)*temp2 - Imatrix*alpha1(:,2) - Fmatrix*alpha2(:,2) - Rmatrix*alpha3(:,2) - Pmatrix.*Imatrix*tau1(:,2) - Pmatrix.*Fmatrix*tau2(:,2) - Pmatrix.*Rmatrix*tau3(:,2);
% r(3) =  rcons(:,3) + delta*rrain(:,3)*rain2 + delta*rtemp(:,3)*temp2 - Imatrix*alpha1(:,3) - Fmatrix*alpha2(:,3) - Rmatrix*alpha3(:,3) - Pmatrix.*Imatrix*tau1(:,3) - Pmatrix.*Fmatrix*tau2(:,3) - Pmatrix.*Rmatrix*tau3(:,3);
% 
% 
% y(1) = x(1) + rcons(:,1) + delta*rrain(:,1)*rain2 + delta*rtemp(:,1)*temp2 - Imatrix*alpha1(:,1) - Fmatrix*alpha2(:,1) - Rmatrix*alpha3(:,1) - Pmatrix.*Imatrix*tau1(:,1) - Pmatrix.*Fmatrix*tau2(:,1) - Pmatrix.*Rmatrix*tau3(:,1);
% y(2) = x(2) + rcons(:,2) + delta*rrain(:,2)*rain2 + delta*rtemp(:,2)*temp2 - Imatrix*alpha1(:,2) - Fmatrix*alpha2(:,2) - Rmatrix*alpha3(:,2) - Pmatrix.*Imatrix*tau1(:,2) - Pmatrix.*Fmatrix*tau2(:,2) - Pmatrix.*Rmatrix*tau3(:,2);
% y(3) = x(3) + rcons(:,3) + delta*rrain(:,3)*rain2 + delta*rtemp(:,3)*temp2 - Imatrix*alpha1(:,3) - Fmatrix*alpha2(:,3) - Rmatrix*alpha3(:,3) - Pmatrix.*Imatrix*tau1(:,3) - Pmatrix.*Fmatrix*tau2(:,3) - Pmatrix.*Rmatrix*tau3(:,3);
% 
