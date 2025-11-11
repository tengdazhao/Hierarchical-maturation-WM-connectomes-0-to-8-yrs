function [indmin mdall notenough]=multicenter_multimodelnlm_global_novalid(Age_gender,Brainsz,data_r,ylabelname,Plotornot,betarealinitial,Site_info)
notenough=[];
Age=Age_gender(:,1)./30./12;
Gender=Age_gender(:,2);
BrainSize=Brainsz./1000;
Age2=Age.*Age;
%fitting data using linear model to obtain the initial beta
if(size(Age_gender,2)>2)
    
Covall=Age_gender(:,3:end);
tbl=table(Age,Gender,BrainSize,Covall,data_r);
tb2=table(Age,Age2,Gender,BrainSize,Covall,data_r);
mdl0line = fitlm(tbl);
mdl0qua = fitlm(tb2);
rpall=size(2,7);
% linear model:mdl1
%fitting data using nonlinear model(fitnlm) with initial beta from fitlm,
%this is because that the nlm is too sensitivy for the initial data and
%easy to fail
beta1 = [mdl0line.Coefficients.Estimate];
modelfun1 = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,2) + b(4)*x(:,3)+b(5)*x(:,4);
mdl1 = fitnlm(tbl,modelfun1,beta1);

% Plot linear model:mdl1
md4plot=mdl1;
tb4plot=tbl;
Data_rremain=data_r-md4plot.Coefficients.Estimate(1)-(md4plot.Coefficients.Estimate(end-2).*Gender)-(md4plot.Coefficients.Estimate(end-1).*BrainSize)-(md4plot.Coefficients.Estimate(end).*Covall);
Data_rremain = Data_rremain+(mean(data_r)-mean(Data_rremain));
[fitresult1, gof1] = fit(Age,Data_rremain,'poly1');
name=[ylabelname 'apoly1Covsp.tif'];
rpall(1,1) = max(md4plot.Coefficients.pValue(2:end-3));
rpall(2,1) = gof1.adjrsquare;
if(Plotornot==1)
    plot4fit_multicenter(Age,Data_rremain,fitresult1,name,ylabelname,rpall(1,1),rpall(2,1),Site_info)
end



% Qua model:mdl2
beta2 = [mdl0qua.Coefficients.Estimate];
modelfun2 = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,1).^2 + b(4)*x(:,2) + b(5)*x(:,3) + b(6)*x(:,4);
mdl2 = fitnlm(tbl,modelfun2,beta2);

% Plot Qua model:mdl2
md4plot=mdl2;
tb4plot=tbl;
Data_rremain=data_r-md4plot.Coefficients.Estimate(1)-(md4plot.Coefficients.Estimate(end-2).*Gender)-(md4plot.Coefficients.Estimate(end-1).*BrainSize)-(md4plot.Coefficients.Estimate(end).*Covall);
Data_rremain = Data_rremain+(mean(data_r)-mean(Data_rremain));
[fitresult2, gof2] = fit(Age,Data_rremain,'poly2');
name=[ylabelname 'bpoly2Covsp.tif'];
rpall(1,2) = max(md4plot.Coefficients.pValue(2:end-3));
rpall(2,2) = gof2.adjrsquare;
if(Plotornot==1)
plot4fit_multicenter(Age,Data_rremain,fitresult2,name,ylabelname,rpall(1,2),rpall(2,2),Site_info)
end

%nonlinear model:
if((rpall(2,1)<0.05)&&(rpall(2,2)<0.05))
    mdl3=mdl1;mdl3=mdl1;mdl5=mdl1;mdl6=mdl1;mdl7=mdl1;
    fitresulte3=fitresult2;fitresul1possion=fitresult2;fitresultlog=fitresult2;fitresultsig1=fitresult2;fitresultsig2=fitresult2;
    rpall(1,3:4)=rpall(1,2);rpall(2,3:4)=rpall(2,2);

else
% exp model:mdl4
% exp model:mdl3
modelnlmreal = @(b,x)b(1) + b(2).*exp(-b(3).*x(:,1))+ b(4)*x(:,2)+ b(5)*x(:,3) +b(6)*x(:,4);
modelfitreal='a.*exp((-b).*x)+c';
Lowerboundary = [-Inf 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
Upperboundary=[Inf 1 Inf];
name=[ylabelname 'cexp3Covsp.tif'];
[mdl3 fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm_multicenter(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial,Site_info);
seor=1;
if(sum(mdl3.Coefficients.pValue==0)>0)
bbnew=~betarealinitial;
    [mdl3new fitresulte3new gofitenew rpalltemp]=Curve_fitting_nlm_multicenter(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,bbnew,Site_info);
[sem seor]=min([mdl3.ModelCriterion.AICc mdl3new.ModelCriterion.AICc]);
end

if(seor==1)
mdl3new=[];[mdl3 fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm_multicenter(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial,Site_info);
else
   mdl3=mdl3new;fitresulte3=fitresulte3new;gofite=gofitenew;rpall(:,3)=rpalltemp;
end
end

else
tbl=table(Age,Gender,BrainSize,data_r);
tb2=table(Age,Age2,Gender,BrainSize,data_r);
mdl0line = fitlm(tbl);
mdl0qua = fitlm(tb2);
rpall=size(2,7);
% linear model:mdl1
%fitting data using nonlinear model(fitnlm) with initial beta from fitlm,
%this is because that the nlm is too sensitivy for the initial data and
%easy to fail
beta1 = [mdl0line.Coefficients.Estimate];
modelfun1 = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,2) + b(4)*x(:,3);
mdl1 = fitnlm(tbl,modelfun1,beta1);

% Plot linear model:mdl1
md4plot=mdl1;
tb4plot=tbl;
Data_rremain=data_r-md4plot.Coefficients.Estimate(1)-(md4plot.Coefficients.Estimate(end-1).*Gender)-(md4plot.Coefficients.Estimate(end).*BrainSize);
Data_rremain = Data_rremain+(mean(data_r)-mean(Data_rremain));
[fitresult1, gof1] = fit(Age,Data_rremain,'poly1');
name=[ylabelname 'apoly1.tif'];
% modelnlmreal = @(b,x)b(1) + b(2).*x(:,1);
% Lowerboundary = [-Inf -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf Inf];
% tblNew=table(Age,Data_rremain);
% modelnlmrealnew = modelnlmreal;
% betarealinitialnew=md4plot.Coefficients.Estimate(1:2);
% mdlnlmnew = fitnlm(tblNew,modelnlmrealnew,betarealinitialnew);
% nobs = mdlnlmnew.NumObservations;
% ssr = max(mdlnlmnew.SST - mdlnlmnew.SSE,0);
% dfr = mdlnlmnew.NumEstimatedCoefficients - 1;
% dfe = nobs - 1 - dfr;
% f = (ssr./dfr) / (mdlnlmnew.SSE/dfe);
% p = fcdf(1./f,dfe,dfr); % upper tail
% Pvalue = max(md4plot.Coefficients.pValue(2:end-2));
% Pvalue = p;
% Rsqaure = gof1.adjrsquare;
rpall(1,1) = max(md4plot.Coefficients.pValue(2:end-2));
rpall(2,1) = gof1.adjrsquare;
if(Plotornot==1)
    plot4fit_multicenter(Age,Data_rremain,fitresult1,name,ylabelname,rpall(1,1),rpall(2,1),Site_info)
end



% Qua model:mdl2
beta2 = [mdl0qua.Coefficients.Estimate];
modelfun2 = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,1).^2 + b(4)*x(:,2) + b(5)*x(:,3);
mdl2 = fitnlm(tbl,modelfun2,beta2);

% Plot Qua model:mdl2
md4plot=mdl2;
tb4plot=tbl;
Data_rremain=data_r-md4plot.Coefficients.Estimate(1)-(md4plot.Coefficients.Estimate(end-1).*Gender)-(md4plot.Coefficients.Estimate(end).*BrainSize);
Data_rremain = Data_rremain+(mean(data_r)-mean(Data_rremain));
[fitresult2, gof2] = fit(Age,Data_rremain,'poly2');
name=[ylabelname 'bpoly2.tif'];
rpall(1,2) = max(md4plot.Coefficients.pValue(2:end-2));
rpall(2,2) = gof2.adjrsquare;
if(Plotornot==1)
plot4fit_multicenter(Age,Data_rremain,fitresult2,name,ylabelname,rpall(1,2),rpall(2,2),Site_info)
end

%nonlinear model:
if((rpall(2,1)<0.05)&&(rpall(2,2)<0.05))
    mdl3=mdl1;mdl3=mdl1;mdl5=mdl1;mdl6=mdl1;mdl7=mdl1;
    fitresulte3=fitresult2;fitresul1possion=fitresult2;fitresultlog=fitresult2;fitresultsig1=fitresult2;fitresultsig2=fitresult2;
    rpall(1,3:4)=rpall(1,2);rpall(2,3:4)=rpall(2,2);

else
% exp model:mdl4
% modelnlmreal = @(b,x)b(1) + b(2).*exp(-b(3).*x(:,1))+ b(4)*x(:,2)+ b(5)*x(:,3);
% modelfitreal='a.*exp((-b).*x)+c';
% Lowerboundary = [-Inf 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf 1 Inf];
% name=[ylabelname 'cexp3.tif'];
% [mdl4o fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);

% % modelnlmreal = @(b,x)b(1) + b(2).*exp(x(:,1))+ b(3)*x(:,2)+ b(4)*x(:,3);
% % modelfitreal='a.*exp(x)+b';
% 
% % log model:mdl3
% modelnlmreal = @(b,x)b(1) + b(2).*log(x(:,1))+ b(3)*x(:,2)+ b(4)*x(:,3);
% modelfitreal='a.*(log(x))+b';
% Lowerboundary = [-Inf -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf Inf];
% name=[ylabelname 'clog3.tif'];
% 
% [mdl3 fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm(mdl1,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);
% 
% seor=1;
% if(sum(mdl3.Coefficients.pValue==0)>0)
% bbnew=~betarealinitial;
%     [mdl3new fitresulte3new gofitenew rpalltemp]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,bbnew);
% [sem seor]=min([mdl3.ModelCriterion.AICc mdl3new.ModelCriterion.AICc]);
% end
% 
% if(seor==1)
% mdl3new=[];[mdl3 fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);
% else
%    mdl3=mdl3new;fitresulte3=fitresulte3new;gofite=gofitenew;rpall(:,3)=rpalltemp;
% end
% 
% exp model:mdl3
modelnlmreal = @(b,x)b(1) + b(2).*exp(-b(3).*x(:,1))+ b(4)*x(:,2)+ b(5)*x(:,3);
modelfitreal='a.*exp((-b).*x)+c';
Lowerboundary = [-Inf 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
Upperboundary=[Inf 1 Inf];
name=[ylabelname 'cexp3.tif'];
[mdl3 fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm_multicenter(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial,Site_info);
seor=1;
if(sum(mdl3.Coefficients.pValue==0)>0)
bbnew=~betarealinitial;
    [mdl3new fitresulte3new gofitenew rpalltemp]=Curve_fitting_nlm_multicenter(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,bbnew,Site_info);
[sem seor]=min([mdl3.ModelCriterion.AICc mdl3new.ModelCriterion.AICc]);
end

if(seor==1)
mdl3new=[];[mdl3 fitresulte3 gofite rpall(:,3)]=Curve_fitting_nlm_multicenter(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial,Site_info);
else
   mdl3=mdl3new;fitresulte3=fitresulte3new;gofite=gofitenew;rpall(:,3)=rpalltemp;
end

%     
% % possion model:mdl4
% modelnlmreal = @(b,x)b(1) + (b(2).*x(:,1)).*exp(-b(3).*x(:,1))+ b(4)*x(:,2)+ b(5)*x(:,3);
% modelfitreal='(a.*x).*exp((-b).*x)+c';
% if(mdl2.Coefficients.Estimate(3)<=0)
% Lowerboundary = [0 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf 1 Inf];
% elseif(mdl2.Coefficients.Estimate(3)>0)
% Lowerboundary = [-Inf 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[0 1 Inf];
% end
% name=[ylabelname 'possion4.tif'];
% [mdl4 fitresul1possion gofitpossion rpall(:,4)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);% real fitting data 
% ttt=0;tta=0;seor=1;
% if(sum(mdl4.Coefficients.pValue==0)>0)
% bbnew=~betarealinitial;
%     [mdl4new fitresul1possionnew gofitpossionnew rpalltemp]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,bbnew);
% [sem seor]=min([mdl4.ModelCriterion.AICc mdl4new.ModelCriterion.AICc]);
% end
% 
% if(seor==1)
% mdl4new=[];[mdl4 fitresul1possion gofitpossion rpall(:,4)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);% real fitting data 
% 
% else
%    mdl4=mdl4new;fitresul1possion=fitresul1possionnew;gofitpossion=gofitpossionnew;rpall(:,4)=rpalltemp;
% end
% while((mdl4.Coefficients.Estimate(3))>1.2)
%     [mdl4 fitresul1possion gofitpossion rpall(:,4)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);
%     ttt=ttt+1;
%     if(ttt>5)
%         break;
%     end
% end


% if(mdl2.Coefficients.Estimate(3)<=0)
% % log model:mdl5
% modelnlmreal = @(b,x)b(1) + b(2).*(log10(x(:,1))./log10(b(3)))+ b(4)*x(:,2)+ b(5)*x(:,3);
% modelfitreal='a.*(log10(x)./log10(b))+c';
% Lowerboundary = [-Inf 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf 1 Inf];
% name=[ylabelname 'elog5.tif'];
% [mdl5 fitresultlog gofitlog rpall(:,5)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot,betarealinitial);
% %mdl5=mdl1;rpall(:,5)=rpall(:,1);
% elseif(mdl2.Coefficients.Estimate(3)>0)
% mdl5=mdl4;fitresultlog=fitresulte3;rpall(:,5)=rpall(:,3);
% end
% % 



% % sigmoid model:mdl6
% modelnlmreal = @(b,x)b(1) + b(2)./(1+exp(-b(3).*x(:,1)))+ b(4)*x(:,2)+ b(5)*x(:,3);
% modelfitreal='a./(1+exp(-b.*x))+c';
% Lowerboundary = [-Inf 0 -Inf];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf 1 Inf];
% name=[ylabelname 'fsigmoid6.tif'];
% [mdl6 fitresultsig1 gofitsig1 rpall(:,6)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot);
% 
% 
% % sigmoid model:mdl7
% modelnlmreal = @(b,x)b(1) + b(2)./(1+exp(-1.*x(:,1)))+ b(3)*x(:,2)+ b(4)*x(:,3);
% modelfitreal='a./(1+exp(-1.*x))+b';
% Lowerboundary = [-Inf -Inf ];%To constrain the lower boundary of beta2 above ZERO.
% Upperboundary=[Inf Inf];
% name=[ylabelname 'gsigmoid7.tif'];
% [mdl7 fitresultsig2 gofitsig2 rpall(:,7)]=Curve_fitting_nlm(mdl2,modelnlmreal,modelfitreal,Lowerboundary,Upperboundary,ylabelname,name,Plotornot);
% end
end
end
[minvalue(1,:) indmin(1,:)]=sort([mdl1.ModelCriterion.BIC mdl3.ModelCriterion.BIC]);
[minvalue(1,:) indmin(1,:)]=sort([mdl1.ModelCriterion.AICc mdl3.ModelCriterion.AICc]);
[indmin(2,1)]=([mdl1.ModelCriterion.AICc]);
[indmin(2,2)]=([mdl3.ModelCriterion.AICc]);
%[indmin(2,3)]=([mdl4.ModelCriterion.AICc]);
[indmin(3,1)]=([mdl1.ModelCriterion.CAIC]);
[indmin(3,2)]=([mdl3.ModelCriterion.CAIC]);
%[indmin(3,3)]=([mdl4.ModelCriterion.CAIC]);
[indmin(4,1)]=([mdl1.ModelCriterion.BIC]);
[indmin(4,2)]=([mdl3.ModelCriterion.BIC]);
%[indmin(4,3)]=([mdl4.ModelCriterion.BIC]);




[indmin(5,:)]=rpall(1,[1 3]);
[indmin(6,:)]=rpall(2,[1 3]);

mdall={mdl1 mdl3;fitresult1 fitresulte3};

if((abs(minvalue(1,1)-minvalue(1,2)))<2)
notenough.ture=1;
notenough.select=min(indmin(1,1:2));
else
notenough.ture=0;
notenough.select=indmin(1,1);  
end

end