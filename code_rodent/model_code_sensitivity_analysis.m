%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB code for sensitivity analysis described in article
%%% The main file for figS5, figS6 and figS7
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
% avg_cul_patch: a column vector of length timePts with average agriculture patch area data
% avg_urb_patch: a column vector of length timePts with average urban patch area data
% temp_dailymax: a column vector of length timePts with daily maximum temperature data

%%%%%%%%%%%%%%%%%%%sensitivity analysis%%%%%%%%%%%%%%%%%%

clear; close all; warning off;
dataset = 'dataset_sensitivity.mat'
%seed(123)
% load(dataset);

% Parameters for fitting the model:
nknots = 10;  % number of knots in the spline fit
h = 1;  % number of smoothing bandwidth (h) 
penalty = -5; % number of penalty weights (penalty) 


%% Sensitivity analysis for substituting mean patch size by agriculture patch area in the three-species dynamic model.
%%%%%%%%%%%%%%%   AA   %%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat,  tau1, tau2, tau3,tau1hat,tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,AAdensity1,AAhat, r2_logAA, r2_AA, r2_GR,  exitflag] = fit_model_AA_cul(dataset, nknots, h, penalty, 1);

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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logAA'); ylabel('logAAhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logAA','logAAhat');

% now get confidence intervals for these results
generateConfidenceIntervalsAA_cul(dataset);

%%%%%%%%%%%%%%%%    RN    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat,  alpha2hat, alpha3hat,tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat,Y,Yhat,RNdensity1,RNhat, r2_logRN,r2_RN,r2_GR, exitflag2] = fit_model_RN_cul(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRN'); ylabel('logRNhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRN','logRNhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRN_cul(dataset);

%%%%%%%%%%%%%%%%    RF    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,RFdensity1,RFhat, r2_logRF,r2_RF,r2_GR, exitflag3] = fit_model_RF_cul(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRF'); ylabel('logRFhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRF','logRFhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRF_cul(dataset);

%% Sensitivity analysis for substituting mean patch size by urban patch area in the three-species dynamic model.
%%%%%%%%%%%%%%%   AA   %%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat,  tau1, tau2, tau3,tau1hat,tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,AAdensity1,AAhat, r2_logAA, r2_AA, r2_GR,  exitflag] = fit_model_AA_urb(dataset, nknots, h, penalty, 1);

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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logAA'); ylabel('logAAhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logAA','logAAhat');

% now get confidence intervals for these results
generateConfidenceIntervalsAA_urb(dataset);

%%%%%%%%%%%%%%%%    RN    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat,  alpha2hat, alpha3hat,tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat,Y,Yhat,RNdensity1,RNhat, r2_logRN,r2_RN,r2_GR, exitflag2] = fit_model_RN_urb(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRN'); ylabel('logRNhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRN','logRNhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRN_urb(dataset);

%%%%%%%%%%%%%%%%    RF    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,RFdensity1,RFhat, r2_logRF,r2_RF,r2_GR, exitflag3] = fit_model_RF_urb(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRF'); ylabel('logRFhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRF','logRFhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRF_urb(dataset);

%% Sensitivity analysis for substituting mean temperature by daily maximum temperature in the three-species dynamic model.
%%%%%%%%%%%%%%%   AA   %%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat,  tau1, tau2, tau3,tau1hat,tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,AAdensity1,AAhat, r2_logAA, r2_AA, r2_GR,  exitflag] = fit_model_AA_tempmax(dataset, nknots, h, penalty, 1);

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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logAA'); ylabel('logAAhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logAA','logAAhat');

% now get confidence intervals for these results
generateConfidenceIntervalsAA_tempmax(dataset);

%%%%%%%%%%%%%%%%    RN    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat,  alpha2hat, alpha3hat,tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat,Y,Yhat,RNdensity1,RNhat, r2_logRN,r2_RN,r2_GR, exitflag2] = fit_model_RN_tempmax(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRN'); ylabel('logRNhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRN','logRNhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRN_tempmax(dataset);

%%%%%%%%%%%%%%%%    RF    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, GRhat, Y,Yhat,RFdensity1,RFhat, r2_logRF,r2_RF,r2_GR, exitflag3] = fit_model_RF_tempmax(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
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

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRF'); ylabel('logRFhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRF','logRFhat');

% now get confidence intervals for these results
generateConfidenceIntervalsRF_tempmax(dataset);