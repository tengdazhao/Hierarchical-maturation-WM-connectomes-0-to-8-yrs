function [corr_stats]=corrstat(data1,data2,cov)
if(~isempty(cov))
N = length(data1);

whichstats={'rsquare','tstat','mse','r'};
% %linear
stats_linear1= regstats(data1,[data2 cov],'linear',whichstats);
beta1 = stats_linear1.tstat.beta;
data_fitted1 = beta1(1)+beta1(2)*data2+stats_linear1.r;


stats_linear2 = regstats(data2,[data1 cov],'linear',whichstats);
beta2= stats_linear2.tstat.beta;
data_fitted2 = beta2(1)+beta2(2)*data1+stats_linear2.r;

% %linear
% stats_linear1= regstats(data1,cov,'linear',whichstats);
% beta1 = stats_linear1.tstat.beta;
% data_fitted1 = beta1(1)+stats_linear1.r;
% 
% 
% stats_linear2 = regstats(data2,cov,'linear',whichstats);
% beta2= stats_linear2.tstat.beta;
% data_fitted2 = beta2(1)+stats_linear2.r;
meandata = mean(data2);
[r p]=partialcorr(data1,data2,cov);
corr_stats = {r,data1,data_fitted2,p,beta2(2),meandata,meandata,beta2(2)./meandata};%4 factor 
else
% corr_stats = {r,data_fitted1,data_fitted2,p};%4 factor
N = length(data1);

whichstats={'rsquare','tstat','mse','r'};
% %linear
stats_linear1= regstats(data1,[data2],'linear',whichstats);
beta1 = stats_linear1.tstat.beta;
data_fitted1 = beta1(1)+beta1(2)*data2+stats_linear1.r;


stats_linear2 = regstats(data2,[data1],'linear',whichstats);
beta2= stats_linear2.tstat.beta;
data_fitted2 = beta2(1)+beta2(2)*data1+stats_linear2.r;
meandata = mean(data2);
% %linear
% stats_linear1= regstats(data1,cov,'linear',whichstats);
% beta1 = stats_linear1.tstat.beta;
% data_fitted1 = beta1(1)+stats_linear1.r;
% 
% 
% stats_linear2 = regstats(data2,cov,'linear',whichstats);
% beta2= stats_linear2.tstat.beta;
% data_fitted2 = beta2(1)+stats_linear2.r;

[r p]=corr(data1,data2);
corr_stats = {r,data1,data_fitted2,p,beta2(2),meandata,beta2(2)./meandata};%4 factor 
end







