function [S] = getSplineSmoother(nknots, m, penalty)

x = (0:m)';

% cubic spline
Xfit(:,1) = ones(length(x),1);
Xfit(:,2) = x;
Xfit(:,3) = x.^2;
Xfit(:,4) = x.^3;

% setting the knots
spacing = length(x)/(nknots + 1);  % evenly spaced knots
for i=1:nknots
    knots(i) = floor(spacing*i);
end

loc = 5;
for i=1:nknots
    Xfit(:,loc) = (max(x-knots(i),0)).^3;
    loc = loc + 1;
end

XTX = Xfit'*Xfit;

% setting the penalty
D = zeros(size(nknots+4));
if penalty ~= -Inf
    for i=1:nknots
        D(i+4, i+4) = 10^penalty;
    end
end    

% SPLINE SMOOTHER- does not include constraints
S = Xfit*inv(Xfit'*Xfit + D)*Xfit';
