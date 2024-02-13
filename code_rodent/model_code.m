%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB code to fit the model described in article
%%% The main file for fig4
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

clear; close all; warning off;
dataset = 'dataset.mat'
%seed(123)

% Parameters for fitting the model:
nknots = 10;  % number of knots in the spline fit
h = 1;  % number of smoothing bandwidth (h) 
penalty = -5; % number of penalty weights (penalty) 

%%%%%%%%%%%%%%%%    AA    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3,tau1hat,tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,AAdensity1,AAhat, r2_logAA, r2_AA, r2_GR,  exitflag] = fit_model_AA(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
figure;
subplot(3,3,1); hold on; ylabel('tau1'); plot(0:(length(tau1)-1), tau1, 'b'); 
subplot(3,3,2); hold on; ylabel('tau2'); plot(0:(length(tau2)-1), tau2, 'b'); 
subplot(3,3,3); hold on; ylabel('tau3'); plot(0:(length(tau3)-1), tau3, 'b'); 
subplot(3,3,4); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(3,3,5); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(3,3,6); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 
subplot(3,3,7); hold on; plot(0:(length(alpha1)-1), alpha1,'b'); ylabel('alpha1'); 
subplot(3,3,8); hold on; plot(0:(length(alpha2)-1), alpha2,'b'); ylabel('alpha2');
subplot(3,3,9); hold on; plot(0:(length(alpha3)-1), alpha3,'b'); ylabel('alpha3'); 

% % get growth rate
% figure;
% subplot(2,1,1); hold on; plot(real(Y0), real(GRhat),'b.');xlabel('GR'); ylabel('GRhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(real(Y0),'LineWidth',1);hold on;plot(real(GRhat),'LineWidth',1,'color','red'); legend('GR','GRhat');

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logAA'); ylabel('logAAhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logAA','logAAhat');

% % get population density
% figure;
% subplot(2,1,1); hold on; plot(real(AAdensity1), real(AAhat),'b.');xlabel('AAdensity'); ylabel('AAdensityhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(real(AAdensity1),'LineWidth',1);hold on;plot(real(AAhat),'LineWidth',1,'color','red'); legend('AAdensity','AAdensityhat');

% now get confidence intervals for these results
generateConfidenceIntervalsAA(dataset);

%%%%%%%%%%%%%%%%    RN    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat,alpha2hat, alpha3hat,  tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat,Y,Yhat,RNdensity1,RNhat, r2_logRN,r2_RN,r2_GR, exitflag2] = fit_model_RN(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file

% plotting the results
figure;
subplot(3,3,1); hold on; ylabel('tau1'); plot(0:(length(tau1)-1), tau1, 'b'); 
subplot(3,3,2); hold on; ylabel('tau2'); plot(0:(length(tau2)-1), tau2, 'b'); 
subplot(3,3,3); hold on; ylabel('tau3'); plot(0:(length(tau3)-1), tau3, 'b'); 
subplot(3,3,4); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(3,3,5); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(3,3,6); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 
subplot(3,3,7); hold on; plot(0:(length(alpha1)-1), alpha1,'b'); ylabel('alpha1'); 
subplot(3,3,8); hold on; plot(0:(length(alpha2)-1), alpha2,'b'); ylabel('alpha2');
subplot(3,3,9); hold on; plot(0:(length(alpha3)-1), alpha3,'b'); ylabel('alpha3'); 

% % get growth rate
% figure;
% subplot(2,1,1); hold on; plot(real(Y0), real(GRhat),'b.');xlabel('GR'); ylabel('GRhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(real(Y0),'LineWidth',1);hold on;plot(real(GRhat),'LineWidth',1,'color','red'); legend('GR','GRhat');

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRN'); ylabel('logRNhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRN','logRNhat');

% % get population density
% figure;
% subplot(2,1,1); hold on; plot(real(RNdensity1), real(RNhat),'b.');xlabel('RNdensity'); ylabel('RNdensityhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(real(RNdensity1),'LineWidth',1);hold on;plot(real(RNhat),'LineWidth',1,'color','red'); legend('RNdensity','RNdensityhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRN(dataset);


%%%%%%%%%%%%%%%%    RF    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,RFdensity1,RFhat, r2_logRF,r2_RF,r2_GR, exitflag3] = fit_model_RF(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
% plotting the results
figure;
subplot(3,3,1); hold on; ylabel('tau1'); plot(0:(length(tau1)-1), tau1, 'b'); 
subplot(3,3,2); hold on; ylabel('tau2'); plot(0:(length(tau2)-1), tau2, 'b'); 
subplot(3,3,3); hold on; ylabel('tau3'); plot(0:(length(tau3)-1), tau3, 'b'); 
subplot(3,3,4); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(3,3,5); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(3,3,6); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 
subplot(3,3,7); hold on; plot(0:(length(alpha1)-1), alpha1,'b'); ylabel('alpha1'); 
subplot(3,3,8); hold on; plot(0:(length(alpha2)-1), alpha2,'b'); ylabel('alpha2');
subplot(3,3,9); hold on; plot(0:(length(alpha3)-1), alpha3,'b'); ylabel('alpha3'); 

% % get growth rate
% figure;
% subplot(2,1,1); hold on; plot(real(Y0), real(GRhat),'b.');xlabel('GR'); ylabel('GRhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(real(Y0),'LineWidth',1);hold on;plot(real(GRhat),'LineWidth',1,'color','red'); legend('GR','GRhat');

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRF'); ylabel('logRFhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRF','logRFhat');

% % get population density
% figure;
% subplot(2,1,1); hold on; plot(real(RFdensity1), real(RFhat),'b.');xlabel('RFdensity'); ylabel('RFdensityhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(real(RFdensity1),'LineWidth',1);hold on;plot(real(RFhat),'LineWidth',1,'color','red'); legend('RFdensity','RFdensityhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRF(dataset);

