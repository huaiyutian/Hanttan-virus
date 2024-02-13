function [r,y] = beddingtonAA (x,rcons,rrain,rtemp,alpha1,tau1,rain2,temp2,delta,Pmatrix,Imatrix)

r(1) =  rcons(:,1) + delta*rrain(:,1)*rain2 + delta*rtemp(:,1)*temp2 - Imatrix*alpha1(:,1)  - Pmatrix.*Imatrix*tau1(:,1); 


y(1) = x(1) + rcons(:,1) + delta*rrain(:,1)*rain2 + delta*rtemp(:,1)*temp2 - Imatrix*alpha1(:,1)  - Pmatrix.*Imatrix*tau1(:,1);
