%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB code to fit the model described in article
%%% The main file for figS4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through smoothing bandwidths (h) and penalty weights (penalty).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supplied dataset contains the following parameters:
% timePts: the number of data points in the dataset
% startTime: the year of the first data point
% n: the number of time intervals that comprise an annual cycle
% m: the maximum duration of time lag
% AAdensity: a column vector of length timePts with population density data of striped field mouse (Apodemus agrarius, AA)
% RNdensity: a column vector of length timePts with population density data of Norway rat (Rattus norvegicus, RN)
% RFdensity: a column vector of length timePts with population density data of buff-breasted rat (R. flavipectus, RF)
% avg_patch_size: a column vector of length timePts with average patch size data
% temp: a column vector of length timePts with temperature data
% rainfall: a column vector of length timePts with rainfall data

%Landscape-ecological consolidation speed affects rodent growth rate.
clear; close all; warning off;
dataset = 'dataset.mat'

% Parameters for fitting the model:
nknots = 10;  % number of knots in the spline fit
h = 1;  % number of smoothing bandwidth (h) 
penalty = -5; % number of penalty weights (penalty) 

%AA
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat,tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,AAdensity1,AAhat, r2_logAA, r2_AA, r2_GR, exitflag] = fit_model_AA(dataset, nknots, h, penalty, 1);

eval(['rcons', '=betalongterm']);
eval(['rrain','=logbetarain']);
eval(['rtemp','=logbetatemp']);
eval(['a1','=alpha1']);
eval(['a2','=alpha2']);
eval(['a3','=alpha3']);
eval(['ttau1','=tau1hat']);
eval(['ttau2','=tau2hat']);
eval(['ttau3','=tau3hat']);

%RN
 [dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat,Y,Yhat,RNdensity1,RNhat,  r2_logRN,r2_RN,r2_GR, exitflag2] = fit_model_RN(dataset, nknots, h, penalty, 1);

a1 = [a1,alpha1];a2 = [a2,alpha2];a3=[a3,alpha3];
ttau1=[ttau1,tau1hat];ttau2=[ttau2,tau2hat];ttau3=[ttau3,tau3hat];
rcons=[rcons,betalongterm];rrain=[rrain,logbetarain];rtemp=[rtemp,logbetatemp];

%RF
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat,Y,Yhat,RFdensity1,RFhat, r2_logRF,r2_RF,r2_GR, exitflag3] = fit_model_RF(dataset, nknots, h, penalty, 1);

a1 = [a1,alpha1];a2 = [a2,alpha2];a3=[a3,alpha3];
ttau1=[ttau1,tau1hat];ttau2=[ttau2,tau2hat];ttau3=[ttau3,tau3hat];
rcons=[rcons,betalongterm];rrain=[rrain,logbetarain];rtemp=[rtemp,logbetatemp];

labels = [{'AAparas'},{'RNparas'},{'RFparas'}];

save('theta.mat','rcons','rrain','rtemp','a1','a2','a3','ttau1','ttau2','ttau3','labels');
clc;clear;
load('theta.mat')
eval(['alpha1','=a1']);
eval(['alpha2','=a2']);
eval(['alpha3','=a3']);
eval(['tau1','=ttau1']);
eval(['tau2','=ttau2']);
eval(['tau3','=ttau3']);
save('theta.mat','rcons','rrain','rtemp','alpha1','alpha2','alpha3','tau1','tau2','tau3','labels');

%%%%%%%%%simulation%%%%%%%%%%
clc;clear;
load('theta.mat')
dataset = 'dataset.mat'
load(dataset);
m=12;
dates(1) = startTime;
for i=2:timePts
    dates(i) = dates(i-1) + 1/n;
end

% using spline inpolatation for AA
nonNaNIndices = find(~isnan(AAdensity));
AAdensity = interp1(nonNaNIndices, AAdensity(nonNaNIndices), 1:length(AAdensity), 'spline');
AAdensity = AAdensity';
% using spline inpolatation for RN
nonNaNIndices = find(~isnan(RNdensity));
RNdensity = interp1(nonNaNIndices, RNdensity(nonNaNIndices), 1:length(RNdensity), 'spline');
RNdensity = RNdensity';
% using spline inpolatation for RF
nonNaNIndices = find(~isnan(RFdensity));
RFdensity = interp1(nonNaNIndices, RFdensity(nonNaNIndices), 1:length(RFdensity), 'spline');
RFdensity = RFdensity';

% % during time points for which there are no cases, replace 0 cases with 1 case (there should not be many of these data points)
locs1 = find(AAdensity <0); AAdensity(locs1)=0;
locs2 = find(RNdensity <0); RNdensity(locs2)=0;
locs3 = find(RFdensity <0); RFdensity(locs3)=0;

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
% delta(:,n) = [];


% generating the vector log(AA(t))
logAA = log(AAdensity);
% replace all infinite values with 0
logAA(isinf(logAA)) = 0;

logRN = log(RNdensity);
logRN(isinf(logRN)) = 0;

logRF = log(RFdensity);
logRF(isinf(logRF)) = 0;

% generating the vector log(AA(t+1))
YAA = zeros(timePts,1);
for i=1:(timePts-1)
	YAA(i,1) = log(AAdensity(i+1,1));
end
YAA(isinf(YAA)) = 0;

% generating the vector log(RN(t+1))
YRN = zeros(timePts,1);
for i=1:(timePts-1)
	YRN(i,1) = log(RNdensity(i+1,1));
end
YRN(isinf(YRN)) = 0;

% generating the vector log(RF(t+1))
YRF = zeros(timePts,1);
for i=1:(timePts-1)
	YRF(i,1) = log(RFdensity(i+1,1));
end
YRF(isinf(YRF)) = 0;

AAdensity1 = zeros(timePts,1);
for i=1:(timePts-1)
	AAdensity1(i,1) = AAdensity(i+1,1);
end

RNdensity1 = zeros(timePts,1);
for i=1:(timePts-1)
	RNdensity1(i,1) = RNdensity(i+1,1);
end

RFdensity1 = zeros(timePts,1);
for i=1:(timePts-1)
	RFdensity1(i,1) = RFdensity(i+1,1);
end

%generating the vector rainfall(t-2)
rain2 = zeros(timePts,1);
for i=1:(timePts-2)
    rain2(i+2,1) = rainfall(i,1);
end

%generating the vector temp(t-2)
temp2 = zeros(timePts,1);
for i=1:(timePts-2)
    temp2(i+2,1) = temp(i,1);
end

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
% Imatrix = Imatrix(:,2:end);

% generating the matrix of previous RNdensity 
Rmatrix = zeros(timePts, m+1);
Rmatrix(:,1) = RNdensity;
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Rmatrix(row,prev+1) = Rmatrix(row-1,prev);
	end
	prev = prev+1;
end
% Rmatrix = Rmatrix(:,2:end);

% generating the matrix of previous RFdensity 
Fmatrix = zeros(timePts, m+1);
Fmatrix(:,1) = RFdensity;
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Fmatrix(row,prev+1) = Fmatrix(row-1,prev);
	end
	prev = prev+1;
end
% Fmatrix = Fmatrix(:,2:end);

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
% Pmatrix = Pmatrix(:,2:end);

% cutting out the first m+1 data points and the last data point for reasons of incompleteness
dates(1:(m)) = [];
rain2([1:(m)],:) = [];
temp2([1:(m)],:) = [];
delta([1:(m)],:) = [];
logAA([1:(m)],:) = [];
logRN([1:(m)],:) = [];
logRF([1:(m)],:) = [];
YAA([1:(m)],:) = [];
YRN([1:(m)],:) = [];
YRF([1:(m)],:) = [];
AAdensity1([1:(m)],:) = [];
RNdensity1([1:(m)],:) = [];
RFdensity1([1:(m)],:) = [];
Imatrix([1:(m)],:) = [];
Rmatrix([1:(m)],:) = [];
Fmatrix([1:(m)],:) = [];
Pmatrix([1:(m)],:) = [];

% resetting the length of the timePts
timePts = timePts-m;

logRodent = [YAA, YRN, YRF];
Rodentdensity1 = [AAdensity1, RNdensity1, RFdensity1];
GR = [YAA-logAA, YRN-logRN, YRF-logRF];

y0GR = beddingtonGR (rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix,Imatrix,Rmatrix,Fmatrix);
% y0logX = beddingtonlogX (logRodent,timePts, y0GR);
% y0X = beddingtonX (Rodentdensity1,timePts, y0GR);

% get growth rate
figure;
subplot(3,1,1);hold on;plot(real(GR(:,1)),'LineWidth',1);hold on;plot(real(y0GR(:,1)),'LineWidth',1,'color','red'); legend('AAGR','AAGRhat');
subplot(3,1,2);hold on;plot(real(GR(:,2)),'LineWidth',1);hold on;plot(real(y0GR(:,2)),'LineWidth',1,'color','red'); legend('RNGR','RNGRhat');
subplot(3,1,3);hold on;plot(real(GR(:,3)),'LineWidth',1);hold on;plot(real(y0GR(:,3)),'LineWidth',1,'color','red'); legend('RFGR','RFGRhat');

% % get log population density
% figure;
% subplot(3,1,1);hold on;plot(real(logRodent(:,1)),'LineWidth',1);hold on;plot(real(y0logX(:,1)),'LineWidth',1,'color','red'); legend('logAA','logAAhat');
% subplot(3,1,2);hold on;plot(real(logRodent(:,2)),'LineWidth',1);hold on;plot(real(y0logX(:,2)),'LineWidth',1,'color','red'); legend('logRN','logRNhat');
% subplot(3,1,3);hold on;plot(real(logRodent(:,3)),'LineWidth',1);hold on;plot(real(y0logX(:,3)),'LineWidth',1,'color','red'); legend('logRF','logRFhat');

% % get population density
% figure;
% subplot(3,1,1);hold on;plot(real(Rodentdensity1(:,1)),'LineWidth',1);hold on;plot(real(y0X(:,1)),'LineWidth',1,'color','red'); legend('AAdensity','AAdensityhat');
% subplot(3,1,2);hold on;plot(real(Rodentdensity1(:,2)),'LineWidth',1);hold on;plot(real(y0X(:,2)),'LineWidth',1,'color','red'); legend('RNdensity','RNdensityhat');
% subplot(3,1,3);hold on;plot(real(Rodentdensity1(:,3)),'LineWidth',1);hold on;plot(real(y0X(:,3)),'LineWidth',1,'color','red'); legend('RFdensity','RFdensityhat');



%% different change rate of patch size
timePts=timePts+m;

% set landscape-ecological consolidation speed
v_change=0.05;

% accelerated land consolidation by 3v
avg_patch_size3=avg_patch_size*(1+3*v_change);

Pmatrix3 = zeros(timePts, m+1);
Pmatrix3(:,1) = avg_patch_size3(1:timePts);
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix3(row,prev+1) = Pmatrix3(row-1,prev);
	end
	prev = prev+1;
end

% accelerated land consolidation by 2v
avg_patch_size2=avg_patch_size*(1+2*v_change);

Pmatrix2 = zeros(timePts, m+1);
Pmatrix2(:,1) = avg_patch_size2(1:timePts);
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix2(row,prev+1) = Pmatrix2(row-1,prev);
	end
	prev = prev+1;
end

% accelerated land consolidation by v
avg_patch_size1=avg_patch_size*(1+1*v_change);

Pmatrix1 = zeros(timePts, m+1);
Pmatrix1(:,1) = avg_patch_size1(1:timePts);
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix1(row,prev+1) = Pmatrix1(row-1,prev);
	end
	prev = prev+1;
end

% slow down land consolidation by v
avg_patch_size1n=avg_patch_size*(1-1*v_change);

Pmatrix1n = zeros(timePts, m+1);
Pmatrix1n(:,1) = avg_patch_size1n(1:timePts);
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix1n(row,prev+1) = Pmatrix1n(row-1,prev);
	end
	prev = prev+1;
end

% slow down land consolidation by 2v
avg_patch_size2n=avg_patch_size*(1-2*v_change);

Pmatrix2n = zeros(timePts, m+1);
Pmatrix2n(:,1) = avg_patch_size2n(1:timePts);
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix2n(row,prev+1) = Pmatrix2n(row-1,prev);
	end
	prev = prev+1;
end

% slow down land consolidation by 3v
avg_patch_size3n=avg_patch_size*(1-3*v_change);
 
Pmatrix3n = zeros(timePts, m+1);
Pmatrix3n(:,1) = avg_patch_size3n(1:timePts);
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix3n(row,prev+1) = Pmatrix3n(row-1,prev);
	end
	prev = prev+1;
end

Pmatrix3([1:(m)],:) = [];
Pmatrix2([1:(m)],:) = [];
Pmatrix1([1:(m)],:) = [];
Pmatrix1n([1:(m)],:) = [];
Pmatrix2n([1:(m)],:) = [];
Pmatrix3n([1:(m)],:) = [];

timePts = timePts-m;

y3GR = beddingtonGR (rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix3,Imatrix,Rmatrix,Fmatrix)
y2GR = beddingtonGR  ( rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix2,Imatrix,Rmatrix,Fmatrix)
y1GR = beddingtonGR  ( rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix1,Imatrix,Rmatrix,Fmatrix)
y1nGR = beddingtonGR  ( rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix1n,Imatrix,Rmatrix,Fmatrix)
y2nGR = beddingtonGR  ( rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix2n,Imatrix,Rmatrix,Fmatrix)
y3nGR = beddingtonGR  ( rcons,rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2,temp2,delta,Pmatrix3n,Imatrix,Rmatrix,Fmatrix)

% y3logX = beddingtonlogX (logRodent,timePts, y3GR)
% y2logX = beddingtonlogX  (logRodent,timePts, y2GR)
% y1logX = beddingtonlogX  (logRodent,timePts, y1GR)
% y1nlogX = beddingtonlogX  (logRodent,timePts, y1nGR)
% y2nlogX = beddingtonlogX  (logRodent,timePts, y2nGR)
% y3nlogX = beddingtonlogX  (logRodent,timePts, y3nGR)
% 
% y3X = beddingtonX (Rodentdensity1,timePts, y3GR)
% y2X = beddingtonX  (Rodentdensity1,timePts, y2GR)
% y1X = beddingtonX  (Rodentdensity1,timePts, y1GR)
% y1nX = beddingtonX  (Rodentdensity1,timePts, y1nGR)
% y2nX = beddingtonX  (Rodentdensity1,timePts, y2nGR)
% y3nX = beddingtonX  (Rodentdensity1,timePts, y3nGR)

yGRfinal=[y3GR,y2GR,y1GR,y0GR,y1nGR,y2nGR,y3nGR];
% ylogXfinal=[y3logX,y2logX,y1logX,y0logX,y1nlogX,y2nlogX,y3nlogX];
% yXfinal=[y3X,y2X,y1X,y0X,y1nX,y2nX,y3nX];


yfinal=yGRfinal;

figure
subplot(3,1,1)
plot(yfinal(:,1))
hold on
for i=1:6
    plot(yfinal(:,3*i+1))
end
hold off
subplot(3,1,2)
plot(yfinal(:,2))
hold on
for i=1:6
    plot(yfinal(:,3*i+2))
end
hold off
subplot(3,1,3)
plot(yfinal(:,3))
hold on
for i=1:6
    plot(yfinal(:,3*i+3))
end
hold off    