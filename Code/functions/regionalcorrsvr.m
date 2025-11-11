function regionalcorrsvr(indnodes,Fitted_Subgroup,Age2,roi_num,Gradeintsall,finalYCorr,typeofstat,result_name,resssname,finalYCorrcov,dirpath);
Rslopecorrall=[];Fitted_ageind=[];tempcolornode=zeros(roi_num,(length(Age2)-1));
Fitted_ageind=1:1:max(Fitted_Subgroup);
Fitted_Subgroupforbin=[ round(length(Age2)/1100) round(length(Age2)/6.95) round(length(Age2)/3.95) round(length(Age2)/2.63) ....
    round(length(Age2)/1.98) round(length(Age2)/1.60) round(length(Age2)/1.30)....
    round(length(Age2)/1.01)];
tempcolornode=[];
tempcolornode(indnodes,Fitted_ageind)=Gradeintsall(:,Fitted_ageind);

finalYCorr2=finalYCorr(find(typeofstat==2));
tempcolornode2=tempcolornode(find(typeofstat==2),:);
finalYCorr1=finalYCorr(find(typeofstat==1));


tempcolornode1=tempcolornode(find(typeofstat==1),:);


tempcolornodeall=tempcolornode(find(typeofstat~=0),:);
finalYCorrall=finalYCorr(find(typeofstat~=0));


if(~isempty(finalYCorrcov))
finalYCorrcov2=finalYCorrcov(find(typeofstat==2));
finalYCorr1cov=finalYCorrcov(find(typeofstat==1));
else
finalYCorr1newcov=[];finalYCorrcov2=[];
end



load([dirpath,'\Her3name.mat']);
finalYCorrSys=[];tempcolornodesys=[];
for i=1:max(AALVonE(:,2));
    indsys=find(AALVonE(:,2)==i);  
    indsyssig=intersect(indnodes,indsys);
    finalYCorrSys{1,i}=finalYCorr(indsyssig);
    tempcolornodesys{1,i}=tempcolornode(indsyssig,:);
end


finalYCorrSysAsso=[finalYCorrSys{2}; finalYCorrSys{3}; finalYCorrSys{7};];
finalYCorrSysPri=[finalYCorrSys{1}; finalYCorrSys{4}; finalYCorrSys{5};finalYCorrSys{6}];

tempcolornodesysAsso=[tempcolornodesys{2}; tempcolornodesys{3}; tempcolornodesys{7}];
tempcolornodesysPri=[tempcolornodesys{1}; tempcolornodesys{4}; tempcolornodesys{5};tempcolornodesys{6}];

[PredictionLOOCVAsso] = SVR_LOOCV(tempcolornodesysAsso, finalYCorrSysAsso', [],'Scale',dirpath);
[PredictionLOOCVPri] = SVR_LOOCV(tempcolornodesysPri, finalYCorrSysPri', [],'Scale',dirpath);
[corr_stats]=corrstat(finalYCorrSysAsso,PredictionLOOCVAsso.Score',[]);
corrplot(finalYCorrSysAsso,corr_stats,'Actual nodal efficiency','Predicted nodal efficiency',[],[],[],[],[],[result_name,'_',resssname,'Predition8yrs.tif'])



%%%first plot
for ii=1:1000
[Prediction_NforderAsso{ii}] = SVR_NFolds(tempcolornodesysAsso, finalYCorrSysAsso',[],2,'Scale',dirpath);
Prediction_NforderAssoCorr(ii,:)=Prediction_NforderAsso{ii}.Corr;

[Prediction_NforderPri{ii}] = SVR_NFolds(tempcolornodesysPri, finalYCorrSysPri', [],2,'Scale',dirpath);
Prediction_NforderPriCorr(ii,:)=Prediction_NforderPri{ii}.Corr;
ii
end
save Prediction_NforderAssoCorr1000ztd.mat Prediction_NforderAssoCorr Prediction_NforderPriCorr

Prediction_NforderAssoCorrarr=Prediction_NforderAssoCorr(:);
histogram(Prediction_NforderAssoCorrarr)
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',45,'LineWidth',3,'FontWeight','bold');
set(gcf,'position', [4         20       1100         990]);
set(gca,'position',[0.15,0.22,0.75,0.75]);
xlabel('Prediction accuracy');
ylabel('Counts');
print(gcf,'-dtiff','-r300',[result_name,'_',resssname,'ASSO_2folders_Rdistrztd','.tif']); close all;



Prediction_NforderPriCorrarr=Prediction_NforderPriCorr(:);

sysname=[];
sysname{1}=['Assoc nodes'];sysname{2}=['Prim nodes'];
meanall=[];


x=[];g=[];

x = [Prediction_NforderAssoCorrarr Prediction_NforderPriCorrarr];
gt1 = repmat({[sysname{1}]},length(Prediction_NforderAssoCorrarr),1);
gt2 = repmat({[sysname{2}]},length(Prediction_NforderPriCorrarr),1);
g=[gt1; gt2];


figure

violinplot(x,g,'GroupOrder',sysname,'ShowData',false,'ShowMean',true,'ShowNotches',true);
set(gca,'XTick',[1:1:length(sysname)],'XTickLabel',sysname);
ylabel('Prediction accuracy');
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',45,'LineWidth',3,'FontWeight','bold');
set(gcf,'position', [4 20  600 1100]);

print(gcf,'-dtiff','-r300',[result_name 'Prediction_Accuarcy_voilinof2types.tif']); close all;



%Weights distibution
PredictionLOOCVAsso=[];Accnew=[];
indexinterp=[1:10:1196];
tempcolornodesysAsso_interp=tempcolornodesysAsso(:,indexinterp);
for ii=1:(size(tempcolornodesysAsso_interp,2));
[PredictionLOOCVAsso{ii}] = SVR_LOOCV(tempcolornodesysAsso_interp(:,1:ii), finalYCorrSysAsso', [],'Scale',dirpath);
Accnew(ii)=PredictionLOOCVAsso{ii}.Corr;
end
save Accnewnewinter10ztd.mat Accnew PredictionLOOCVAsso


Fitted_Subgroupforbin=[ round(length(Age2)/1100) round(length(Age2)*(1/16)) round(length(Age2)*(2/16))...
    round(length(Age2)*(3/16)) round(length(Age2)*(4/16))...
    round(length(Age2)*(5/16)) round(length(Age2)*(6/16))...
    round(length(Age2)*(7/16)) round(length(Age2)*(8/16))...
    round(length(Age2)*(9/16)) round(length(Age2)*(10/16))...
    round(length(Age2)*(11/16)) round(length(Age2)*(12/16))...
    round(length(Age2)*(13/16)) round(length(Age2)*(14/16))...
    round(length(Age2)*(15/16)) 1180];


Fitted_Subgroupforbintemp=Fitted_Subgroupforbin;
Fitted_Subgroupforbinnew=round(Fitted_Subgroupforbintemp./10)+1;
w_Brain = Accnew;
w_Brainabs = abs(w_Brain);

w_Brainabs(w_Brain < 0) = 0;
data = w_Brainabs;
Agenew = [1:1:length(w_Brainabs)];
binEdges = [Fitted_Subgroupforbinnew]; 


binIndices = discretize(Agenew, binEdges);
binnum = max(binIndices);


binMeans = nan(1, binnum); 
binmeanind = nan(1, binnum);
for i = binnum:-1:1 
    binIdx = find(binIndices == i);
    if ~isempty(binIdx)
        binMeans(i) = mean(data(binIdx));
        binmeanind(i) = max(binIdx);
    end
end


figure;
barh(binMeans); 
title('Cumulative distribution');
xlabel('Prediction accuracy'); 
ylabel('Ages'); 
maxBinMean = max(binMeans);
set(gca, 'XTick', [0.1 0.3 0.5 0.7], ... 
    'XTickLabel', [0.1 0.3 0.5 0.7], ... 
    'YTick', 1:2:binnum, ... 
    'YTickLabel', [0 1 2 3 4 5 6 7 8]); 
set(gca, 'YDir', 'reverse')
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'fontsize', 35, 'LineWidth', 3, 'FontWeight', 'bold');
set(gcf, 'position', [4, 20, 1350, 900]);
set(gca, 'position', [0.18, 0.28, 0.65, 0.55]);

print('-dpng', '-r300', 'NgE_FDR_SVR_Ageat8ysACCumulation_distribution_ASSO.tif');
close all
% 
% data=w_Brainabs;
% Agenew=[1:1:length(w_Brainabs)]
% Fitted_Subgroupforbin=Fitted_Subgroupforbin
% binEdges = [Fitted_Subgroupforbin];
% binIndices = discretize(Agenew, binEdges);
% binnum=max(binIndices);
% 
% binMeans=[];binmeanind=[];
% for i = 1:binnum
%     binMeans(i) = mean(data(find(binIndices==i)));
%     binmeanind(i)=max(find(binIndices==i));
% end
% figure
% bar(binMeans);
% title('Prediction weights distribution');
% Age2(Fitted_Subgroupforbin)
% xlabel('Ages');
% ylabel('Prediction weights');
% set(gca,'XTick',[1:1:binnum],'XTickLabel',[1 2 3 4 5 6 8]);
% 
% set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',55,'LineWidth',3,'FontWeight','bold');
% set(gcf,'position', [4         20       1450         1000]);
% 
% print(gcf,'-dtiff','-r300',[result_name,'_',resssname,'weightsalldistribution_ASSO','.tif']); close all;
% close all
% corrname=[result_name,'AssonodesPredict','LOOCV',resssname,'.tif'];
% corrstats=corrstat(finalYCorrSysAsso,PredictionLOOCVAsso.Score',[]);corrplot(finalYCorrSysAsso,corrstats,xtitlelabel, ytitlelabel,[],[],[],[],[],corrname);





randnn=10000;
Accnew_rand=zeros(randnn,2);
PredictionLOOCVAsso=[];PredictionLOOCVPri=[];
for jj=1:randnn
finalYCorrSysAsso_random=finalYCorrSysAsso(randperm(length(finalYCorrSysAsso)));
finalYCorrSysPri_random=finalYCorrSysPri(randperm(length(finalYCorrSysPri)));
% for ii=1:(length(tempcolornodesysAsso));
% [PredictionLOOCVAsso{ii}] = SVR_LOOCV(tempcolornodesysAsso(:,1:ii), finalYCorrSysAsso_random', [],'Scale',dirpath);
% Accnew_rand(ii,jj)=PredictionLOOCVAsso{ii}.Corr;
% end
[PredictionLOOCVAsso] = SVR_LOOCV(tempcolornodesysAsso, finalYCorrSysAsso_random', [],'Scale',dirpath);
Accnew_rand(jj,1)=PredictionLOOCVAsso.Corr;
[PredictionLOOCVPri] = SVR_LOOCV(tempcolornodesysPri, finalYCorrSysPri_random', [],'Scale',dirpath);
Accnew_rand(jj,2)=PredictionLOOCVPri.Corr;
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
jj
end
save Accnew_rand10000ztd.mat Accnew_rand



