clc;clear;

xdata1 = xlsread('Hu 2023.xlsx');
xdata1 = xdata1(1:516,:);
AAdensity = xdata1(:,3);
RNdensity = xdata1(:,4);
RFdensity = xdata1(:,5);
avg_patch_size = xdata1(:,6);
rainfall = xdata1(:,7);
temp = xdata1(:,8);
m = 12;
startTime = 1980;
n = 12;
timePts = 516;

save('dataset.mat', 'AAdensity','RFdensity','RNdensity','avg_patch_size','temp','rainfall','m','n','startTime','timePts');

xdata1 = xlsread('Hu 2023.xlsx');
xdata1 = xdata1(1:516,:);
AAdensity = xdata1(:,3);
RNdensity = xdata1(:,4);
RFdensity = xdata1(:,5);
avg_patch_size = xdata1(:,6);
rainfall = xdata1(:,7);
temp = xdata1(:,8);
m = 12;
startTime = 1980;
n = 12;
timePts = 516;
cul_area = xdata1(:,12);
cul_num = xdata1(:,13);
avg_cul_patch = cul_area./cul_num;
urb_area = xdata1(:,14);
urb_num = xdata1(:,15);
avg_urb_patch = cul_area./cul_num;
temp_dailymax = xdata1(:,16);
reshaped_data = reshape(temp_dailymax(49:516), [12, 39]); 
monthly_means = mean(reshaped_data, 2);  
temp_dailymax(1:48) = repmat(monthly_means, [1, 4]);

save('dataset_sensitivity.mat', 'AAdensity','RFdensity','RNdensity','avg_patch_size',...
    'temp','rainfall','m','n','startTime','timePts','avg_cul_patch','avg_urb_patch','temp_dailymax');
									
