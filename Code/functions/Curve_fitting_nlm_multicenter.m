function [mdlnlm fitresult gofit rp]=Curve_fitting_nlm_multicenter(initialmd,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betainitial,Site_info)
tbl=initialmd.Variables;%Using fit to constrain the lower boundary of beta2 above ZERO (Or the model will be failed).
if(size(tbl,2)==5)
Data_rremain=tbl.data_r-initialmd.Coefficients.Estimate(1)-(initialmd.Coefficients.Estimate(end-2).*tbl.Gender)-(initialmd.Coefficients.Estimate(end-1).*tbl.BrainSize)-(initialmd.Coefficients.Estimate(end).*tbl.Covall);
Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
elseif(size(tbl,2)==6)
Data_rremain=tbl.data_r-initialmd.Coefficients.Estimate(1)-(initialmd.Coefficients.Estimate(end-3).*tbl.Gender)-(initialmd.Coefficients.Estimate(end-2).*tbl.BrainSize)-(initialmd.Coefficients.Estimate(end-1).*tbl.Covall1)-(initialmd.Coefficients.Estimate(end).*tbl.Covall2);
Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
else
Data_rremain=tbl.data_r-initialmd.Coefficients.Estimate(1)-(initialmd.Coefficients.Estimate(end-1).*tbl.Gender)-(initialmd.Coefficients.Estimate(end).*tbl.BrainSize);
Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
end
ftmodel = fittype(modelfitreal, 'independent', 'x', 'dependent', 'y' );opts = fitoptions(ftmodel);

opts.Algorithm = 'Trust-Region';
%opts.Algorithm = 'Levenberg-Marquardt';
opts.DiffMaxChange = 0.2;
opts.Display = 'Off';
opts.MaxFunEvals = 6000;
opts.MaxIter = 10000;
opts.Lower = Lowerboundary;%To constrain the lower boundary of beta2 above ZERO.
opts.Upper = Upperboundary;
% opts.Lower(end)=initialmd.Coefficients.Estimate(1)-5;
% opts.Upper(end)=initialmd.Coefficients.Estimate(1)+5;
%[ns numberofx]=size(Lowerboundary);
% ForfitbetaEstimate=[initialmd.Coefficients.Estimate(2:numberofx); initialmd.Coefficients.Estimate(1)];
% opts.StartPoint = [ForfitbetaEstimate];
[fitresult_initial gofit_initial] = fit(tbl.Age,Data_rremain,ftmodel,opts);
%tbl.data_r=Data_rremain;
%mdlnlm = fitnlm(tbl,modelnlmreal,betarealinitial);% real fitting data
%  Pvalue = 0;
%  Rsqaure = gofit_initial.adjrsquare;
%  plot4fit(tbl.Age,Data_rremain,fitresult_initial,name,ylabelname,Pvalue,Rsqaure);
betafit=coeffvalues(fitresult_initial);
%obtain the inital beta for covariates
if(size(tbl,2)==5)
betacov=[initialmd.Coefficients.Estimate(end-2) initialmd.Coefficients.Estimate(end-1) initialmd.Coefficients.Estimate(end)];
betamain=[betafit(end) betafit(1:end-1)];
elseif(size(tbl,2)==6)
betacov=[initialmd.Coefficients.Estimate(end-3:end)'];
betamain=[betafit(end) betafit(1:end-1)];
else
betacov=[initialmd.Coefficients.Estimate(end-1) initialmd.Coefficients.Estimate(end)];
betamain=[betafit(end) betafit(1:end-1)];
end
%betarealinitial = [betafit(3) betafit(1) betafit(2) betagender betabs];%fitting data using fitnlm with initial beta from fit(with constrain of beta2)
if(sum(betainitial)==1)
betarealinitial = [betamain betacov];%fitting data using fitnlm with initial beta from fit(with constrain of beta2)
else
betarealinitial = zeros(length([betamain betacov]),1);;
end
mdlnlm = fitnlm(tbl,modelnlmreal,betarealinitial);% real fitting data 
% Plot possion model:mdl3
md4plot=mdlnlm;
if(size(tbl,2)==5)
Data_rremain=tbl.data_r-initialmd.Coefficients.Estimate(1)-(initialmd.Coefficients.Estimate(end-2).*tbl.Gender)-(initialmd.Coefficients.Estimate(end-1).*tbl.BrainSize)-(initialmd.Coefficients.Estimate(end).*tbl.Covall);
Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
elseif(size(tbl,2)==6)
Data_rremain=tbl.data_r-initialmd.Coefficients.Estimate(1)-(initialmd.Coefficients.Estimate(end-3).*tbl.Gender)-(initialmd.Coefficients.Estimate(end-2).*tbl.BrainSize)-(initialmd.Coefficients.Estimate(end-1).*tbl.Covall1)-(initialmd.Coefficients.Estimate(end).*tbl.Covall2);
Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
else
Data_rremain=tbl.data_r-md4plot.Coefficients.Estimate(1)-(md4plot.Coefficients.Estimate(end-1).*tbl.Gender)-(md4plot.Coefficients.Estimate(end).*tbl.BrainSize);
Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
end
Age=tbl.Age;


opts.Algorithm = 'Trust-Region';
%opts.Algorithm = 'Levenberg-Marquardt';
opts.DiffMaxChange = 0.2;
opts.Display = 'Off';
opts.MaxFunEvals = 6000;
opts.MaxIter = 1000;
opts.Lower = Lowerboundary;%To constrain the lower boundary of beta2 above ZERO.
[ns numberofx]=size(Lowerboundary);
opts.Upper = Upperboundary;
ForfitbetaEstimate=[md4plot.Coefficients.Estimate(2:numberofx); md4plot.Coefficients.Estimate(1)];

% ForfitbetaEstimateNew=[mdlnlmnew.Coefficients.Estimate(2:numberofx); mdlnlmnew.Coefficients.Estimate(1)];
% Forbound=[ForfitbetaEstimate ForfitbetaEstimateNew ];
opts.Lower = [ForfitbetaEstimate(1:end-1)'-0.0001 -Inf];%To constrain the lower boundary of beta2 above ZERO.
opts.Upper = [ForfitbetaEstimate(1:end-1)'+0.0001 Inf];

opts.StartPoint = [ForfitbetaEstimate];
[fitresult gofit] = fit(tbl.Age,Data_rremain,ftmodel,opts);

% tblNew=table(Age,Data_rremain);
% charmodel=char(modelnlmreal);
% charmodel=str2func(charmodel(1:end-24));
% modelnlmrealnew = charmodel;
% betarealinitialnew=[ForfitbetaEstimate(end); ForfitbetaEstimate(1:end-1)];
% mdlnlmnew = fitnlm(tblNew,modelnlmrealnew,betarealinitialnew);
% 
% nobs = mdlnlmnew.NumObservations;
% ssr = max(mdlnlmnew.SST - mdlnlmnew.SSE,0);
% dfr = mdlnlmnew.NumEstimatedCoefficients - 1;
% dfe = nobs - 1 - dfr;
% f = (ssr./dfr) / (mdlnlmnew.SSE/dfe);
% p = fcdf(1./f,dfe,dfr); % upper tail
% Pvalue = max(md4plot.Coefficients.pValue(2:end-2));

%[pt]= coefTest(md4plot,[0 1 0 0 0;0 0 1 0 0;]);%F test
if(size(tbl,2)==5)
Pvalue =  max(md4plot.Coefficients.pValue(2:end-3));
elseif(size(tbl,2)==6)
Pvalue =  max(md4plot.Coefficients.pValue(2:end-4));
else
Pvalue =  max(md4plot.Coefficients.pValue(2:end-2));
end
Rsqaure = gofit.adjrsquare;
rp(1,1)=Pvalue;
rp(2,1)=Rsqaure;
if(Plotornot==1)
plot4fit_multicenter(tbl.Age,Data_rremain,fitresult,name,ylabelname,rp(1,1),rp(2,1),Site_info);
end