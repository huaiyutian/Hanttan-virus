function [dates, alpha_val, gamma_val, tau1hat, betaseas, betalongterm, Y, Yhat, r2,  exitflag, cv, exitflagCV] = fit_model(dataset, nknots, h, penalty, saveBool);

load(dataset);

% generating the dates of the time series
dates(1) = startTime;
for i=2:timePts
    dates(i) = dates(i-1) + 1/n;
end

% during time points for which there are no cases, replace 0 cases with 1 case (there should not be many of these data points)
locs = find(AAdensity == 0);
AAdensity(locs) = 1;

% generating the delta matrix
delta = zeros(timePts, n);
for row=1:timePts
    seas = mod(row, n);
    if seas == 0
        seas = n;
    end
    delta(row,seas) = 1;
end
% taking off the last column to make the number of columns in the delta matrix only be of length n-1
delta(:,n) = [];

% generating the vector log(AA(t))
logAA = log(AAdensity);
    
% generating the vector log(AA(t+1))
Y = zeros(timePts,1);
for i=1:(timePts-1)
	Y(i,1) = log(AAdensity(i+1,1));
end

%generating the vector rainfall(t-2)
rain = zeros(timePts,1);
for i=1:(timePts-2)
    rain(i+2,1) = rainfall(i,1);
end

%generating the vector temp(t-2)
temp2 = zeros(timePts,1);
for i=1:(timePts-2)
    temp2(i+2,1) = temp2(i,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% generating the matrix of previous AAdensity 
Imatrix = zeros(timePts, m+1);
Imatrix(:,1) = AAdensity;

prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Imatrix(row,prev+1) = Imatrix(row-1,prev);
	end
	prev = prev+1;
end
% generating the matrix of previous patch effect
Pmatrix = zeros(timePts, m+1);
Pmatrix(:,1) = avg_patch_size;

prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix(row,prev+1) = Pmatrix(row-1,prev);
	end
	prev = prev+1;
end

% cutting out the first m+1 data points and the last data point for reasons of incompleteness
dates(timePts) = [];
dates(1:(m+1)) = [];
rain([timePts],:) = [];
rain([1:(m+1)],:) = [];
temp2([timePts],:) = [];
temp2([1:(m+1)],:) = [];
delta([timePts], :) = [];
delta([1:(m+1)],:) = [];
logAA([timePts], :) = [];
logAA([1:(m+1)],:) = [];
Y([timePts],:) = [];
Y([1:(m+1)],:) = [];
Imatrix([timePts], :) = [];
Imatrix([1:(m+1)],:) = [];
Pmatrix([timePts], :) = [];
Pmatrix([1:(m+1)],:) = [];

% resetting the length of the timePts
timePts = timePts-(m+2);

% generating the long-term beta smoother
Wblt = getweightmatrix(h, timePts); 

% generating an identity matrix of the same dimensions as the long-term beta smoother
Iblt = eye(size(Wblt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have all the necessary vectors and matrices that we'll keep on re-using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepping for the first weighted least-squares regression:

% getting Z
for i=1:timePts
    Z(i,:) = -Pmatrix(i,:).*Imatrix(i,:);
end
% getting delta_rain
delta_rain=zeros(size(delta));
for i=1:timePts
    delta_rain(i,:) = delta(i,:).*rain(i);
end

% getting delta_temp
delta_temp=zeros(size(delta));
for i=1:timePts
    delta_temp(i,:) = delta(i,:).*temp2(i);
end

%getting delta_temp*rain
rain_temp=rain.*temp2;
delta_raintemp=zeros(size(delta));
for i=1:timePts
    delta_raintemp(i,:) = delta(i,:).*rain_temp(i);
end

% getting X
% X = [delta, logAA, Z];
% X = [delta_rain, logAA, Z];
% X = [delta_temp, logAA, Z];
X = [delta_raintemp, logAA, Z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST WEIGHTED-LEAST-SQUARES REGRESSION
theta = inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y;
logbetaseas = [theta([1:n-1],1); 0];
alpha_val= theta(n,1); gamma_val = 1;
tau1 = theta([n+1:(n+m+1)]);
% tau1 = theta([n+1,1]);
% tau2 = theta([n+2,1]);
% tau3 = theta([n+3,1]);
% kappa = theta([n+1:(n+m+1)]);
[tau1hat, exitflag] = getsplinefit(tau1, nknots, m, penalty);
% due to setting gamma_val at 1, the spline-fit kappa may be producing a negative number of susceptibles. getinboundskappa therefore adjusts kappa to make S/N always between 0 and 1.
% [tau1hat, inbounds] = getinboundskappa(Imatrix, tau1hat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBSEQUENT WEIGHTED LEAST-SQUARES REGRESSIONS: The backfitting algorithm

backfit = 1;
% while exitflag ~= 0 & inbounds ~= -1
   
for k =1:20    
    %for i=1:timePts
     for i=1:size(Pmatrix, 1)
        Z(i,:) = -Pmatrix(i,:).*Imatrix(i,:);
    end
    
    X = [delta, logAA, Z*tau1hat, Z];
    theta = inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y;
    logbetaseas = [theta([1:n-1],1); 0];
    alpha_val = theta(n,1); 
    gamma_val = theta(n+1,1);
    gamma_phi = theta([(n+2):(m+n+2)],1);
    phi = gamma_phi./gamma_val;
    
    tau1 = tau1hat + phi;
    %[tau1hat, exitflag] = getsplinefit(tau1, nknots, m, penalty);
%     
%     if(backfit > 20), break, end;  % don't force S/N to be positive anymore on the last backfit. 20 backfits usually are more than sufficient for convergence of kappahat
    
    backfit = backfit+1;
    
%     [tau1hat, inbounds] = getinboundskappa(Imatrix, tau1hat, N);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOW PROCEED WITH FINAL FIT   
if max(tau1hat) == 0
    inbounds = -1;
end

% Y = f1hat + f2hat, where f1hat is the parametric part of the equation and f2hat is the nonparametric part of the equation
X = [delta, logAA]; f1hat = X*theta(1:n,1); 

betaseas = exp(logbetaseas);

% f1hat is the parametric component of the fit
% f2hat is the non-parametric component of the fit

f2hat = Wblt*(Y-f1hat); % local constant fit
betalongterm = exp(f2hat); 

% entire fit.
Yhat = f1hat+f2hat; resids = Y-Yhat;

% R^2
SSR = sum((Yhat-mean(Y)).^2);
SSTO = sum((Y-mean(Y)).^2);
r2 = SSR/SSTO;

% cross-validated weight matrix
Wbltcv = getweightmatrixcv(h, timePts); 
f2hatcv = Wbltcv*(Y-f1hat);

% cross-validated kappa fit
[tau1hatCV, exitflagCV] = getsplinefitcv(tau1, nknots, m, penalty);

% if max(tau1hatCV) == 0
%     inboundsCV = -1;
% end

% putting the cv parts together
Xcv = [delta, logAA]; f1hatcv = Xcv*theta(1:n,1); 
Yhatcv = f1hatcv+f2hatcv; cv = sum((Y-Yhatcv).^2);

if saveBool
    filename = strcat('Opt_Results_', dataset);
    save(filename, 'nknots', 'm', 'n', 'penalty', 'h', 'dates', 'Wblt', 'Iblt', 'X', 'Y', 'Yhat', 'alpha_val', 'gamma_val', 'tau1', 'tau1hat', 'logbetaseas', 'f2hat');
end