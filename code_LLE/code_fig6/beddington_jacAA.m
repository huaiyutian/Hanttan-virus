function y = beddington_jacAA (x,rcons,rrain,rtemp,alpha1,tau1,rain2,temp2,delta,Pmatrix,Imatrix)


% Jac for AA
for j=1:13
    y(1,j) = -alpha1(j,1) - Pmatrix(j)*tau1(j,1);
end
y(1,2) = 1/exp(x(1)) + y(1,2) ;




