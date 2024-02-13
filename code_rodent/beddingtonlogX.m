function y = beddingtonlogX (x, timePts,GR)
y = zeros(timePts, 3);

y(1,1) = x(1,1);
y(1,2) = x(1,2);
y(1,3) = x(1,3);

for i=2:timePts
    y(i,1) = y(i-1,1) + GR(i,1);
    y(i,2) = y(i-1,2) + GR(i,2);
    y(i,3) = y(i-1,3) + GR(i,3);
end 
