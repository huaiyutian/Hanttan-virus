function y = beddingtonX(x, timePts, GR)

y = zeros(timePts, 3);

y(3,1)=x(3,1);
for i = 4:timePts
    y(i,1)=y(i-1,1)*exp(GR(i,1));
end

y(5,2)=x(5,2);
for i = 6:timePts
    y(i,2)=y(i-1,2)*exp(GR(i,2));
end

y(16,3)=x(16,3);
for i = 17:timePts
    y(i,3)=y(i-1,3)*exp(GR(i,3));
end

