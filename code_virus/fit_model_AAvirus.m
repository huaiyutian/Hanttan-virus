function [dates, betap, betaphat, logbetarain,logbetatemp, betalongterm, Y, Yhat,Imatrix1,Ihat, r2,r2_virus1] = fit_model_AAvirus(dataset, nknots, h, saveBool);
load(dataset);
m=12;
rng(123);

% generating the dates of the time series
dates(1) = startTime;
for i=2:timePts
    dates(i) = dates(i-1) + 1/n;
end

% using spline inpolatation for perAAvirus(1,49:468)
nonNaNIndices = find(~isnan(perAAvirus));
perAAvirus = interp1(nonNaNIndices, perAAvirus(nonNaNIndices), 1:length(perAAvirus), 'spline');
perAAvirus = perAAvirus';
noise = 0.01;
% add noise to the original data
perAAvirus = real(perAAvirus + noise);
% during time points for which there are no cases, replace 0 cases with 1 case (there should not be many of these data points)
locs1 = find(perAAvirus >=1); perAAvirus(locs1)=1;
locs2 = find(perAAvirus <0); perAAvirus(locs2)=0;

perAAvirus = round(perAAvirus, 4); 

S=1-perAAvirus;

% generating the delta matrix
delta = zeros(timePts, n);
for row=1:timePts
    seas = mod(row, n);
    if seas == 0
        seas = n;
    end
    delta(row,seas) = 1;
end
% % taking off the last column to make the number of columns in the delta matrix only be of length n-1
delta(:,n) = [];

% generating the vector log(perAAvirus(t))
logAAv = log(perAAvirus);
% replace all infinite values with 0
logAAv(isinf(logAAv)) = 0;

%generating the vector rainfall(t-2)
rain = zeros(timePts,1);
for i=1:(timePts-2)
    rain(i+2,1) = rainfall(i,1);
end

%generating the vector temp(t-2)
temp2 = zeros(timePts,1);
for i=1:(timePts-2)
    temp2(i+2,1) = temp(i,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% generating the matrix of previous AAdensity 
Imatrix = zeros(timePts, 1);
Imatrix(:,1) = perAAvirus;

Imatrix1 = zeros(timePts,1);
for i=1:(timePts-1)
	Imatrix1(i,1) = Imatrix(i+1,1);
end

logImatrix=log(Imatrix);
Smatrix=log(1-Imatrix);
logImatrix(isinf(logImatrix)) = 0;
% generate random noise that follows a standard normal distribution
noise = -0.01* randn(size(logImatrix));
% add noise to the original data
logImatrix = real(logImatrix + noise);
locs3 = find(logImatrix >0); logImatrix(locs3)=-0.01;

Smatrix(isinf(Smatrix)) = 0;
locs4 = find(Smatrix >0); Smatrix(locs4)=-0.01;
locs5 = find(abs(Smatrix) <0.001); Smatrix(locs5)=-0.01;

% generating the vector log(perAAvirus(t+1))
Y = zeros(timePts,1);
for i=1:(timePts-1)
	Y(i,1) = logImatrix(i+1,1);
end
Y(isinf(Y)) = 0;


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
dates(1:(m)) = [];
rain([1:(m)],:) = [];
temp2([1:(m)],:) = [];
delta([1:(m)],:) = [];
logAAv([1:(m)],:) = [];
Y([1:(m)],:) = [];
Imatrix([1:(m)],:) = [];
Imatrix1([1:(m)],:) = [];
logImatrix([1:(m)],:) = [];
Smatrix([1:(m)],:) = [];
Pmatrix([1:(m)],:) = [];
avg_patch_size([1:(m)],:) = [];
S(1:(m)) = [];


% resetting the length of the timePts
timePts = timePts-m;

% generating the long-term beta smoother
Wblt = getweightmatrix(h, timePts); 

% generating an identity matrix of the same dimensions as the long-term beta smoother
Iblt = eye(size(Wblt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have all the necessary vectors and matrices that we'll keep on re-using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

X = real([delta_rain, delta_temp, Pmatrix, logImatrix, Smatrix]);
Y0=Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST WEIGHTED-LEAST-SQUARES REGRESSION
theta = real(inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y0);
logbetarain = [theta([1:n-1],1); 0];
logbetatemp = [theta([n:2*n-2],1); 0];
betap = theta([2*n-1:2*n+m-1],1);
alphap = theta([2*n+m],1);
gamma = theta([2*n+m+1],1);

% Fitting using csaps
x = [1, 2, 3,4,5,6,7,8,9,10,11,12,13];
p1 = csaps(x, betap', 0.8);  % The third parameter is the smoothness parameter, which can be adjusted to control the smoothness of the fitting
betaphat = fnval(p1,x);
betaphat = betaphat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBSEQUENT WEIGHTED LEAST-SQUARES REGRESSIONS: The backfitting algorithm
   
for k =1:50    
    X = [delta_rain, delta_temp,logImatrix, Smatrix];
    Y1 = Y0-Pmatrix*betaphat;
    theta = real(inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y1);
    logbetarain = [theta([1:n-1],1); 0];
    logbetatemp = [theta([n:2*n-2],1); 0];
    alphap = theta([2*n-1],1);
    gamma = theta([2*n],1);
    
    rho = Y1-X*theta;
    Z = Pmatrix;
    
    phi = inv(Z'*(Iblt-Wblt)*Z)*Z'*(Iblt-Wblt)*rho;
    betap = betaphat+phi;
    
    % Fitting using csaps
    x = [1, 2, 3,4,5,6,7,8,9,10,11,12,13];
    p1 = csaps(x, betap', 0.8);  % The third parameter is the smoothness parameter, which can be adjusted to control the smoothness of the fitting
    betaphat = fnval(p1,x);
    betaphat = betaphat';
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOW PROCEED WITH FINAL FIT 
betaphat = real(betaphat);

% Y = f1hat + f2hat, where f1hat is the parametric part of the equation and f2hat is the nonparametric part of the equation

X = [delta_rain, delta_temp,logImatrix, Smatrix];
Y1 = real(Y0-Pmatrix*betaphat);
theta = real(inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y1);
f1hat = X*theta; 

% f1hat is the parametric component of the fit
% f2hat is the non-parametric component of the fit

f2hat = Wblt*(Y1-f1hat); % local constant fit
betalongterm = f2hat; 

% entire fit.
% Yhat = logAA+f1hat+f2hat+U*tau1hat+V*tau2hat+W*tau3hat; 
Yhat = f1hat+f2hat+Pmatrix*betaphat; 
resids = Y-Yhat;
Ihat = exp(Yhat);

% R^2
SSR = sum((Yhat-mean(Y)).^2);
SSTO = sum((Y-mean(Y)).^2);
r2 = SSR/SSTO;

SSR1 = sum((Ihat-mean(Imatrix1)).^2);
SSTO1 = sum((Imatrix1-mean(Imatrix1)).^2);
r2_virus1 = SSR1/SSTO1;

if saveBool
    filename = strcat('Opt_Results_AAvirus_', dataset);
    save(filename, 'nknots', 'm', 'n', 'h', 'dates', 'Wblt', 'Iblt', 'X', 'Y', 'Yhat', 'betap', 'betaphat','logbetarain','logbetatemp','betalongterm', 'f2hat');
end