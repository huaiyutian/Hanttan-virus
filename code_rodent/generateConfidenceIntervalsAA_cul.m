function void = generateConfidenceIntervalsAA_cul(dataset)

filename = strcat('Opt_Results_AA_cul_', dataset)
load(filename);

tau1_errors = tau1-tau1hat;
tau2_errors = tau2-tau2hat;
tau3_errors = tau3-tau3hat;
alpha1_errors = alpha1-alpha1hat;
alpha2_errors = alpha2-alpha2hat;
alpha3_errors = alpha3-alpha3hat;

S = getSplineSmoother(nknots, m, penalty);

warning off
for sample = 1:1000
    for i=1:(m+1)
        rand_values = randperm(m+1);
        tau1_new(i,1) = tau1hat(i) + tau1_errors(rand_values(1));
        tau2_new(i,1) = tau2hat(i) + tau2_errors(rand_values(1));
        tau3_new(i,1) = tau3hat(i) + tau3_errors(rand_values(1));
        alpha1_new(i,1) = alpha1hat(i) + alpha1_errors(rand_values(1));
        alpha2_new(i,1) = alpha2hat(i) + alpha2_errors(rand_values(1));
        alpha3_new(i,1) = alpha3hat(i) + alpha3_errors(rand_values(1));
    end
    [tau1hat_new, exitflag1] = getsplinefit(tau1_new, nknots, m, penalty);
    [tau2hat_new, exitflag2] = getsplinefit2(tau2_new, nknots, m, penalty);
    [tau3hat_new, exitflag3] = getsplinefit3(tau3_new, nknots, m, penalty);
    [alpha1hat_new, exitflag4] = getsplinefit4(alpha1_new, nknots, m, penalty);
    [alpha2hat_new, exitflag5] = getsplinefit5(alpha2_new, nknots, m, penalty);
    [alpha3hat_new, exitflag6] = getsplinefit6(alpha3_new, nknots, m, penalty);
    exitflag1_list(sample) = exitflag1;
    exitflag2_list(sample) = exitflag2;
    exitflag3_list(sample) = exitflag3;
    exitflag4_list(sample) = exitflag4;
    exitflag5_list(sample) = exitflag5;
    exitflag6_list(sample) = exitflag6;
    tau1_new_list(:,sample) = tau1_new;
    tau2_new_list(:,sample) = tau2_new;
    tau3_new_list(:,sample) = tau3_new;
    alpha1_new_list(:,sample) = alpha1_new;
    alpha2_new_list(:,sample) = alpha2_new;
    alpha3_new_list(:,sample) = alpha3_new;
    tau1hat_new_list(:,sample) = tau1hat_new;
    tau2hat_new_list(:,sample) = tau2hat_new;
    tau3hat_new_list(:,sample) = tau3hat_new;
    alpha1hat_new_list(:,sample) = alpha1hat_new;
    alpha2hat_new_list(:,sample) = alpha2hat_new;
    alpha3hat_new_list(:,sample) = alpha3hat_new;
    varianceStar1(sample)= getBootstrapVariance(tau1_new, nknots, m);
    varianceStar2(sample)= getBootstrapVariance(tau2_new, nknots, m);
    varianceStar3(sample)= getBootstrapVariance(tau3_new, nknots, m);
    varianceStar4(sample)= getBootstrapVariance(alpha1_new, nknots, m);
    varianceStar5(sample)= getBootstrapVariance(alpha2_new, nknots, m);
    varianceStar6(sample)= getBootstrapVariance(alpha3_new, nknots, m);
    vStar1_g(sample) = (tau1hat_new-tau1hat)'*inv(S*S'*varianceStar1(sample))*(tau1hat_new-tau1hat);
    vStar2_g(sample) = (tau2hat_new-tau2hat)'*inv(S*S'*varianceStar2(sample))*(tau2hat_new-tau2hat);
    vStar3_g(sample) = (tau3hat_new-tau3hat)'*inv(S*S'*varianceStar3(sample))*(tau3hat_new-tau3hat);
    vStar4_g(sample) = (alpha1hat_new-alpha1hat)'*inv(S*S'*varianceStar4(sample))*(alpha1hat_new-alpha1hat);
    vStar5_g(sample) = (alpha2hat_new-alpha2hat)'*inv(S*S'*varianceStar5(sample))*(alpha2hat_new-alpha2hat);
    vStar6_g(sample) = (alpha3hat_new-alpha3hat)'*inv(S*S'*varianceStar6(sample))*(alpha3hat_new-alpha3hat);
end

[sortedVStar1_g, sortedIndexList1] = sort(vStar1_g);
[sortedVStar2_g, sortedIndexList2] = sort(vStar2_g);
[sortedVStar3_g, sortedIndexList3] = sort(vStar3_g);
[sortedVStar4_g, sortedIndexList4] = sort(vStar4_g);
[sortedVStar5_g, sortedIndexList5] = sort(vStar5_g);
[sortedVStar6_g, sortedIndexList6] = sort(vStar6_g);
% confidence interval = 95%, 80%, etc.
confInt = 0.95;
removeLocs1 = (floor(confInt*length(sortedVStar1_g))+1):length(sortedVStar1_g);
sortedVStar1_g(removeLocs1) = [];
sortedIndexList1(removeLocs1) = [];
removeLocs2 = (floor(confInt*length(sortedVStar2_g))+1):length(sortedVStar2_g);
sortedVStar2_g(removeLocs2) = [];
sortedIndexList2(removeLocs2) = [];
removeLocs3 = (floor(confInt*length(sortedVStar3_g))+1):length(sortedVStar3_g);
sortedVStar3_g(removeLocs3) = [];
sortedIndexList3(removeLocs3) = [];
removeLocs4 = (floor(confInt*length(sortedVStar4_g))+1):length(sortedVStar4_g);
sortedVStar4_g(removeLocs4) = [];
sortedIndexList4(removeLocs4) = [];
removeLocs5 = (floor(confInt*length(sortedVStar5_g))+1):length(sortedVStar5_g);
sortedVStar5_g(removeLocs5) = [];
sortedIndexList5(removeLocs5) = [];
removeLocs6 = (floor(confInt*length(sortedVStar6_g))+1):length(sortedVStar6_g);
sortedVStar6_g(removeLocs6) = [];
sortedIndexList6(removeLocs6) = [];

randPermList1 = randperm(length(sortedVStar1_g));
shuffledList1 = sortedVStar1_g(randPermList1);
shuffledIndexList1 = sortedIndexList1(randPermList1);
randPermList2 = randperm(length(sortedVStar2_g));
shuffledList2 = sortedVStar2_g(randPermList2);
shuffledIndexList2 = sortedIndexList2(randPermList2);
randPermList3 = randperm(length(sortedVStar3_g));
shuffledList3 = sortedVStar3_g(randPermList3);
shuffledIndexList3 = sortedIndexList3(randPermList3);
randPermList4 = randperm(length(sortedVStar4_g));
shuffledList4 = sortedVStar4_g(randPermList4);
shuffledIndexList4 = sortedIndexList4(randPermList4);
randPermList5 = randperm(length(sortedVStar5_g));
shuffledList5 = sortedVStar5_g(randPermList5);
shuffledIndexList5 = sortedIndexList5(randPermList5);
randPermList6 = randperm(length(sortedVStar6_g));
shuffledList6 = sortedVStar6_g(randPermList6);
shuffledIndexList6 = sortedIndexList6(randPermList6);

optLoc1 = shuffledIndexList1(1:10);
optLoc2 = shuffledIndexList2(1:10);
optLoc3 = shuffledIndexList3(1:10);
optLoc4 = shuffledIndexList4(1:10);
optLoc5 = shuffledIndexList5(1:10);
optLoc6 = shuffledIndexList6(1:10);

figure; subplot(3,3,4); hold off; plot(0:m, tau1hat_new_list(:,optLoc1), 'k'); hold on; plot(0:m, tau1hat, 'b'); line(0:m, 0);ylabel('tau1hat');
subplot(3,3,5); hold on; plot(0:m, tau2hat_new_list(:,optLoc2), 'k'); hold on; plot(0:m, tau2hat, 'b'); line(0:m, 0);ylabel('tau2hat');
subplot(3,3,6); hold on; plot(0:m, tau3hat_new_list(:,optLoc3), 'k'); hold on; plot(0:m, tau3hat, 'b'); line(0:m, 0);ylabel('tau3hat');
subplot(3,3,7); hold on; plot(0:m, alpha1hat_new_list(:,optLoc4), 'k'); hold on; plot(0:m, alpha1hat, 'b'); line(0:m, 0);ylabel('alpha1hat');
subplot(3,3,8); hold on; plot(0:m, alpha2hat_new_list(:,optLoc5), 'k'); hold on; plot(0:m, alpha2hat, 'b'); line(0:m, 0);ylabel('alpha2hat');
subplot(3,3,9); hold on; plot(0:m, alpha3hat_new_list(:,optLoc6), 'k'); hold on; plot(0:m, alpha3hat, 'b'); line(0:m, 0);ylabel('alpha3hat');

R1 = inv(X'*(Iblt-Wblt)*X);
dferr = length(Y) - trace(2*R1-R1*R1');
MSE = sum((Y-Yhat).^2)/dferr;

var_cov_matrix_est = R1*MSE;
conf_intervals_param_est = 2*sqrt(diag(var_cov_matrix_est));

subplot(3,3,2); hold on;
plot(logbetarain(1:(n-1))+conf_intervals_param_est(1:(n-1)), 'k.'); plot(logbetarain(1:(n-1))-conf_intervals_param_est(1:(n-1)), 'k.');
plot(logbetarain, 'b'); hold on; plot(logbetarain, 'b.');ylabel('betarain');
subplot(3,3,3); hold on;
plot(logbetatemp(1:(n-1))+conf_intervals_param_est(1:(n-1)), 'k.'); plot(logbetatemp(1:(n-1))-conf_intervals_param_est(1:(n-1)), 'k.');
plot(logbetatemp, 'b'); hold on; plot(logbetatemp, 'b.');ylabel('betatemp');

R2 = Wblt*(Iblt-X*inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt));
dferr = length(Y) - trace(2*R2-R2*R2');
MSE = sum((Y-Yhat).^2)/dferr;

cov_f2hat_est = R2*R2'*MSE;
conf_bands_beta_est = 2*sqrt(diag(cov_f2hat_est));

subplot(3,3,1); hold on;
plot(dates, f2hat + conf_bands_beta_est, 'k');
plot(dates, f2hat - conf_bands_beta_est, 'k');
plot(dates, f2hat, 'b');xlabel('time');ylabel('logAAdensity');
