clc;clear;
func = @beddington;
funcjac = @beddington_jac;

load('dataset.mat');
load('theta.mat');
%%
rng(123);
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
locs1 = find(AAdensity <0); AAdensity(locs1)=0;
locs2 = find(RNdensity <0); RNdensity(locs2)=0;
locs3 = find(RFdensity <0); RFdensity(locs3)=0;

%%
% generating the delta matrix
delta = zeros(timePts, n);
for row=1:timePts
    seas = mod(row, n);
    if seas == 0
        seas = n;
    end
    delta(row,seas) = 1;
end

% generating the vector log(AA(t))
logAA = log(AAdensity);
% replace all infinite values with 0
logAA(isinf(logAA)) = 0;
% generating the vector log(RN(t))
logRN = log(RNdensity);
% replace all infinite values with 0
logRN(isinf(logRN)) = 0;
% generating the vector log(RF(t))
logRF = log(RFdensity);
% replace all infinite values with 0
logRF(isinf(logRF)) = 0;

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

rain2([1:m],:) = [];
temp2([1:m],:) = [];
delta([1:m],:) = [];
logAA([1:m],:) = [];
logRN([1:m],:) = [];
logRF([1:m],:) = [];
Imatrix([1:m],:) = [];
Rmatrix([1:m],:) = [];
Fmatrix([1:m],:) = [];
Pmatrix([1:m],:) = [];
AAdensity([1:m],:) = [];
RNdensity([1:m],:) = [];
RFdensity([1:m],:) = [];
timePts = timePts-m;

%% Local Lyapunov expoent
lyap = [];
for j = 2:498
    Imatrix0 = Imatrix(1,:);
    Rmatrix0 = Rmatrix(1,:);
    Fmatrix0 = Fmatrix(1,:);
    x0 = [0,-4.5,-4.3];
    for i = 1:j-1
        t=i+1;
        [r,x] = feval(func,x0, rcons(t-1,:),rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2(t-1,:),temp2(t-1,:),delta(t-1,:),Pmatrix(t-1,:),Imatrix0,Rmatrix0,Fmatrix0);
        AA_t = Imatrix0(1,1)*(1+r(1));
        RN_t = Rmatrix0(1,1)*(1+r(2));
        RF_t = Fmatrix0(1,1)*(1+r(3));
        Imatrix1=Imatrix0;Rmatrix1=Rmatrix0;Fmatrix1=Fmatrix0;
        for i = 2:m
            Imatrix0(1,i) = Imatrix1(1,i-1);
            Rmatrix0(1,i) = Rmatrix1(1,i-1);
            Fmatrix0(1,i) = Fmatrix1(1,i-1);
        end
        Imatrix0(1,1) = AA_t;
        Rmatrix0(1,1) = RN_t ;
        Fmatrix0(1,1) = RF_t;
    end
    Ta=6;
    s=0;
    for t = j+1:(j+6+1)
        [r,x]= feval(func,x,rcons(t-1,:),rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2(t-1,:),temp2(t-1,:),delta(t-1,:),Pmatrix(t-1,:),Imatrix0,Rmatrix0,Fmatrix0);
        AA_t = Imatrix0(1,1)*(1+r(1));
        RN_t = Rmatrix0(1,1)*(1+r(2));
        RF_t = Fmatrix0(1,1)*(1+r(3));
        jac = feval(funcjac,x,rcons(t-1,:),rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2(t-1,:),temp2(t-1,:),delta(t-1,:),Pmatrix(t-1,:),Imatrix0,Rmatrix0,Fmatrix0);
        jac(isinf(jac)) = 0;
        Imatrix1=Imatrix0;Rmatrix1=Rmatrix0;Fmatrix1=Fmatrix0;
        for i = 2:m
            Imatrix0(1,i) = Imatrix1(1,i-1);
            Rmatrix0(1,i) = Rmatrix1(1,i-1);
            Fmatrix0(1,i) = Fmatrix1(1,i-1);
        end
        Imatrix0(1,1) = AA_t;
        Rmatrix0(1,1) = RN_t;
        Fmatrix0(1,1) = RF_t;
        s=s+log(norm(jac,inf));
    end   
    LLE = s/Ta;
    lyap = [lyap; LLE];
end

s=0;
for j = 2:504
    x = [0,-4.5,-4.3];
    Imatrix0 = Imatrix(1,:);
    Rmatrix0 = Rmatrix(1,:);
    Fmatrix0 = Fmatrix(1,:);
    [r,x] = feval(func,x,rcons(t-1,:),rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2(t-1,:),temp2(t-1,:),delta(t-1,:),Pmatrix(t-1,:),Imatrix0,Rmatrix0,Fmatrix0);
    AA_t = Imatrix0(1,1)*(1+r(1));
    RN_t = Rmatrix0(1,1)*(1+r(2));
    RF_t = Fmatrix0(1,1)*(1+r(3));
    jac = feval(funcjac, x, rcons(t-1,:),rrain,rtemp,alpha1,alpha2,alpha3,tau1,tau2,tau3,rain2(t-1,:),temp2(t-1,:),delta(t-1,:),Pmatrix(t-1,:),Imatrix0,Rmatrix0,Fmatrix0);
    Imatrix1=Imatrix0;Rmatrix1=Rmatrix0;Fmatrix1=Fmatrix0;
    for i = 2:m
        Imatrix0(1,i) = Imatrix1(1,i-1);
        Rmatrix0(1,i) = Rmatrix1(1,i-1);
        Fmatrix0(1,i) = Fmatrix1(1,i-1);
    end
    Imatrix0(1,1) = AA_t;
    Rmatrix0(1,1) = RN_t ;
    Fmatrix0(1,1) = RF_t;
    s=s+log(norm(jac,inf));
end
GLE = s/504;
    

figure;
plot(lyap(:,1),'k-','Linewidth', 1);
hold on;
plot([0 500],[GLE GLE],'r');
xlabel('Time','FontName','times','FontSize',16)
ylabel('LLE','FontName','times','FontSize',16)
hold off

