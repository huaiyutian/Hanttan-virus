%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB code to fit the model described in article
%%% The main file for fig5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supplied dataset contains the following parameters:
% timePts: the number of data points in the dataset
% startTime: the year of the first data point
% n: the number of time intervals that comprise an annual cycle
% m: the maximum duration of time lag
% AAdensity: a column vector of length timePts with population density data of striped field mouse (Apodemus agrarius, AA)
% perAAvirus: a column vector of length timePts with Hantaan virus prevalence among striped field mice (AA)
% avg_patch_size: a column vector of length timePts with average patch size data
% temp: a column vector of length timePts with temperature data
% rainfall: a column vector of length timePts with rainfall data

clear; close all; warning off;
dataset = 'dataset_virus.mat'

% Parameters for fitting the model:
nknots = 10;  % number of knots in the spline fit
h = 1;  % number of smoothing bandwidth (h) 

%%%%%%%%%%%%%%%%    AAvirus    %%%%%%%%%%%%%%%%%%%

[dates, betap, betaphat, logbetarain,logbetatemp, betalongterm, Y, Yhat,Imatrix1,Ihat, r2,r2_virus1] = fit_model_AAvirus(dataset, nknots, h, 1);
     
figure;
subplot(2,2,1); hold on; ylabel('logbatap'); plot(betap,'b'); 
subplot(2,2,2); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(2,2,3); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(2,2,4); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 

figure;
plot(0:(length(betap)-1), betap, 'b'); ylabel('patch effect on virus'); xlabel('lags'); 

% figure;
% subplot(2,1,1); hold on; plot(Y, Yhat,'b.');xlabel('logAAvirus'); ylabel('logAAvirushat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
% subplot(2,1,2);hold on;plot(Y,'LineWidth',1);hold on;plot(Yhat,'LineWidth',1,'color','red'); legend('logAAvirus','logAAvirushat');

figure;
subplot(2,1,1); hold on; plot(Imatrix1, Ihat,'b.');xlabel('AAvirus'); ylabel('AAvirushat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(Imatrix1,'LineWidth',1);hold on;plot(Ihat,'LineWidth',1,'color','red'); legend('AAvirus','AAvirushat');

% now get confidence intervals for these results
generateConfidenceIntervalsAAvirus(dataset);



