function [variance] = getBootstrapVariance(tau1, nknots, m)

tau1hat = getsplinefit(tau1, nknots, m, -Inf); % little/no smoothing penalty- according to p. 3 of Ruppert and Carroll

Sinf = getSplineSmoother(nknots, m, -Inf);

variance = (tau1hat-tau1)'*(tau1hat-tau1)/(m-2*trace(Sinf)+trace(Sinf*Sinf'));  % p. 62 Hastie + Tibshirani- same as p. 3 of R + C