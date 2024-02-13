function y = beddington_jac (x,rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix,Imatrix,Rmatrix,Fmatrix)

for j=1:13
    y(1,j) = -alpha1(j,1) - Pmatrix(j)*tau1(j,1);
    y(2,j) = -alpha1(j,2) - Pmatrix(j)*tau1(j,2);
    y(3,j) = -alpha1(j,3) - Pmatrix(j)*tau1(j,3);
end
y(1,2) = 1/exp(x(1)) + y(1,2) ;


for j=1:13
    y(1,13+j) = -alpha2(j,1) - Pmatrix(j)*tau2(j,1);
    y(2,13+j) = -alpha2(j,2) - Pmatrix(j)*tau2(j,2);
    y(3,13+j) = -alpha2(j,3) - Pmatrix(j)*tau2(j,3);
end
y(2,15) = 1/exp(x(2)) + y(2,15) ;


for j=1:13
    y(1,26+j) = -alpha3(j,1) - Pmatrix(j)*tau3(j,1);
    y(2,26+j) = -alpha3(j,2) - Pmatrix(j)*tau3(j,2);
    y(3,26+j) = -alpha3(j,3) - Pmatrix(j)*tau3(j,3);
end
y(3,28) = 1/exp(x(3)) + y(3,28);

% for j=1:13
%     y(1,j) = -alpha1(j,1) - Pmatrix(j)*tau1(j,1);
%     y(2,j) = -alpha1(j,2) - Pmatrix(j)*tau1(j,2);
%     y(3,j) = -alpha1(j,3) - Pmatrix(j)*tau1(j,3);
% end
% y(1,2) = 1/exp(x(1)) + y(1,2) ;
% 
% 
% for j=1:13
%     y(1,13+j) = -alpha2(j,1) - Pmatrix(j)*tau2(j,1);
%     y(2,13+j) = -alpha2(j,2) - Pmatrix(j)*tau2(j,2);
%     y(3,13+j) = -alpha2(j,3) - Pmatrix(j)*tau2(j,3);
% end
% y(2,15) = 1/exp(x(2)) + y(2,15) ;
% 
% 
% for j=1:13
%     y(1,26+j) = -alpha3(j,1) - Pmatrix(j)*tau3(j,1);
%     y(2,26+j) = -alpha3(j,2) - Pmatrix(j)*tau3(j,2);
%     y(3,26+j) = -alpha3(j,3) - Pmatrix(j)*tau3(j,3);
% end
% y(3,28) = 1/exp(x(3)) + y(3,28);




