function void = generateConfidenceIntervalsAAvirus(dataset)

filename = strcat('Opt_Results_AAvirus_', dataset)
load(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIDENCE INTERVALS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

betap_errors = betap-betaphat;

penalty = -5;% number of penalty weights (penalty) 

S = getSplineSmoother(nknots, m, penalty);

warning off
for sample = 1:1000
    for i=1:(m+1)
        rand_values = randperm(m+1);
        betap_new(i,1) = betaphat(i) + betap_errors(rand_values(1));
    end
    % Fitting using csaps
    x = [1, 2, 3,4,5,6,7,8,9,10,11,12,13];
    p1 = csaps(x, betap_new', 0.8);  % The third parameter is the smoothness parameter, which can be adjusted to control the smoothness of the fitting
    betaphat_new = fnval(p1,x);
    betaphat_new = betaphat_new';
    betap_new_list(:,sample) = betap_new;
    betaphat_new_list(:,sample) = betaphat_new;
    varianceStar1(sample)= getBootstrapVariance(betap_new, nknots, m);
    vStar1_g(sample) = (betaphat_new-betaphat)'*inv(S*S'*varianceStar1(sample))*(betaphat_new-betaphat);
end

[sortedVStar1_g, sortedIndexList1] = sort(vStar1_g);
% confidence interval = 95%, 80%, etc.
confInt = 0.95;
removeLocs1 = (floor(confInt*length(sortedVStar1_g))+1):length(sortedVStar1_g);
sortedVStar1_g(removeLocs1) = [];
sortedIndexList1(removeLocs1) = [];

randPermList1 = randperm(length(sortedVStar1_g));
shuffledList1 = sortedVStar1_g(randPermList1);
shuffledIndexList1 = sortedIndexList1(randPermList1);

optLoc1 = shuffledIndexList1(1:10);


figure; hold on; plot(0:m, betaphat_new_list(:,optLoc1), 'k'); hold on; plot(0:m, betaphat, 'b'); line(0:m, 0);ylabel('betaphat');

