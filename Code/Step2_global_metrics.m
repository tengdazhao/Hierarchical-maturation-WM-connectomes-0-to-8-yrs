%Global metrics


% load graph theoritical measurements
load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\DegreeCentrality\DegreeCentrality.mat']);
load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\SmallWorld\SmallWorld.mat']);
load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\NetworkEfficiency\NetworkEfficiency.mat']);
load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);


[roi_num]=size(Matrixallpar500FN{1,1},2);
indall_select=[1:length(Brainsz)]';
plotornot=1;

cd(Pathnewplot)
Sp=[];AvgLength=[];AvgnormLength=[];NodalLength=[];NodalStrst=[];NodalStrlon=[];stpall=[];Matrixfl=[];Matrixf2=[];AvgnormLength=[];



%sp
for i = 1:length(Matrixallpar500FN)
    %Matrix = Fiber_Count_mult_by_FA{i,1};
    Matrix = Matrixallpar500FN{i,1};
    MatrixB = Matrix;
    MatrixB(find(MatrixB>0))=1;
    Sp(i,1) = length(find(MatrixB>0))/length(Matrix)/(length(Matrix)-1);
end

Age_gender_sp=[Age_gender Sp];

[indminEg mdallEg notEg]=multimodelnlm_global(Age_gender_sp,Brainsz,Eg,'Csp Global Eff',plotornot,indall_select,1)
[indminEloc mdallEloc notEloc]=multimodelnlm_global(Age_gender_sp,Brainsz,Eloc,'Csp Local_E',plotornot,indall_select,1)
[indminGamma mdallGamma notGamma]=multimodelnlm_global(Age_gender_sp,Brainsz,Gamma,'Csp Gamma',plotornot,indall_select,1)
[indminLambda mdallLambda notLambda]=multimodelnlm_global(Age_gender_sp,Brainsz,Lambda,'Csp Lambda',plotornot,indall_select,1)
[indminSigma mdallSigma notSigma]=multimodelnlm_global(Age_gender_sp,Brainsz,Sigma,'Csp Sigma',plotornot,indall_select,1)

%robustness
allrob=[];
allrob=dir([Pathnew,'Forgithub\Data\Gretna_results\Par500FN\Robustness\*.mat']);
TargetLargeComSizeall=[];TargetLpall=[];RandomLargeComSizeall=[];RandomLpall=[];
LarComSize=[];Lp=[];
for i =1:size(allrob,1);
load([allrob(i).folder,'\',allrob(i).name]);
TargetLargeComSizeall(:,i) =[LarComSize.targetattack,zeros(1,roi_num-length(LarComSize.targetattack))];
TargetLpall(:,i) =[Lp.targetattack,zeros(1,roi_num-length(Lp.targetattack))];
RandomLargeComSizeall(:,i) =[LarComSize.randomattack,zeros(1,roi_num-length(LarComSize.randomattack))];
RandomLpall(:,i) =[Lp.randomattack,zeros(1,roi_num-length(Lp.randomattack))];
end

TargetEffall=1./TargetLpall;
RandomEffall=1./RandomLpall;
TargetEffall(find(TargetEffall==Inf))=0;
RandomEffall(find(RandomEffall==Inf))=0;

AUC_Eff_target =sum(TargetEffall./500,1);
AUC_Eff_random =sum(RandomEffall./500,1);

AUC_LCC_target =sum(TargetLargeComSizeall./500,1);
AUC_LCC_random =sum(RandomLargeComSizeall./500,1);

[indmintarget mdalltarget notEtarget]=multimodelnlm_global(Age_gender_sp,Brainsz,AUC_Eff_target','Csp AUC Eff target',plotornot,indall_select,1)
[indminrandom mdallrandom notrandom]=multimodelnlm_global(Age_gender_sp,Brainsz,AUC_Eff_random','Csp AUC Eff random',plotornot,indall_select,1)


%growth rates
h=plot(mdallEg{2,indminEg(1,1)},Age,mdallEg{1,indminEg(1,1)}.Fitted, 'predobs');
hxvalue(:,1)=h(2,1).XData;hyvalue(:,1)=h(2,1).YData;hyhvalue(:,1)=h(3,1).YData;hylvalue(:,1)=h(4,1).YData;
h=plot(mdallEloc{2,indminEloc(1,1)},Age,mdallEloc{1,indminEloc(1,1)}.Fitted, 'predobs');hold on;
hxvalue(:,2)=h(2,1).XData;hyvalue(:,2)=h(2,1).YData;hyhvalue(:,2)=h(3,1).YData;hylvalue(:,2)=h(4,1).YData;
h=plot(mdallGamma{2,indminGamma(1,1)},Age,mdallGamma{1,indminGamma(1,1)}.Fitted, 'predobs');hold on;
hxvalue(:,3)=h(2,1).XData;hyvalue(:,3)=h(2,1).YData;hyhvalue(:,3)=h(3,1).YData;hylvalue(:,3)=h(4,1).YData;
h=plot(mdallSigma{2,indminSigma(1,1)},Age,mdallSigma{1,indminSigma(1,1)}.Fitted, 'predobs');hold on;
hxvalue(:,4)=h(2,1).XData;hyvalue(:,4)=h(2,1).YData;hyhvalue(:,4)=h(3,1).YData;hylvalue(:,4)=h(4,1).YData;
h=plot(mdalltarget{2,indmintarget(1,1)},Age,mdalltarget{1,indmintarget(1,1)}.Fitted, 'predobs');hold on;
hxvalue(:,5)=h(2,1).XData;hyvalue(:,5)=h(2,1).YData;hyhvalue(:,5)=h(3,1).YData;hylvalue(:,5)=h(4,1).YData;
close all
Gradeintsall=[];Gradeintshall=[];Gradeintslall=[];
for ii=1:5
Dx=diff(hxvalue(:,ii));
Dy=diff(hyvalue(:,ii));
Dxy=Dy./Dx;
Gradeintsall(:,ii)=Dxy;
Dx=diff(hxvalue(:,ii));
Dyh=diff(hyhvalue(:,ii));
Dxyh=Dyh./Dx;
Gradeintshall(:,ii)=Dxyh;
Dx=diff(hxvalue(:,ii));
Dyl=diff(hylvalue(:,ii));
Dxyl=Dyl./Dx;
Gradeintslall(:,ii)=Dxyl;
ii
end
figure;
yData=Gradeintsall(:);
xData=Age;
axes('LineWidth',5,'FontWeight','bold','FontSize',48);
axis([min(xData)-(max(xData)-min(xData))/10 max(xData)+(max(xData)-min(xData))/12 min(yData)-(max(yData)-min(yData))/15 max(yData)+(max(yData)-min(yData))/15]);
hold on;
hnewh=plot(hxvalue(2:end,:),Gradeintshall);
hold on;
hnewl=plot(hxvalue(2:end,:),Gradeintslall);
hold on;
hnew=plot(hxvalue(2:end,:),Gradeintsall);

set(hnewh(2),'LineWidth',0.5,'LineStyle','-','color',[255/255 255/255 255/255]);
set(hnewh(1),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    254/255]);
set(hnewh(5),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
set(hnewh(3),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
set(hnewh(4),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);

set(hnewl(2),'LineWidth',0.5,'LineStyle','-','color',[255/255 255/255 255/255]);
set(hnewl(1),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
set(hnewl(5),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
set(hnewl(3),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
set(hnewl(4),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
% 
fill([hnewh(2).XData,fliplr(hnewh(2).XData)],[hnewh(2).YData,fliplr(hnewl(2).YData)],[227/255 150/255 120/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色
fill([hnewh(1).XData,fliplr(hnewh(1).XData)],[hnewh(1).YData,fliplr(hnewl(1).YData)],[254/255    210/255    150/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色
fill([hnewh(5).XData,fliplr(hnewh(5).XData)],[hnewh(5).YData,fliplr(hnewl(5).YData)],[120/255    226/255    200/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色
fill([hnewh(3).XData,fliplr(hnewh(3).XData)],[hnewh(3).YData,fliplr(hnewl(3).YData)],[255/255     169/255    160/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色
fill([hnewh(4).XData,fliplr(hnewh(4).XData)],[hnewh(4).YData,fliplr(hnewl(4).YData)],[255/255     169/255    160/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色

set(hnew(2),'LineWidth',2.5,'LineStyle','-','color',[227/255 68/255 10/255]);
set(hnew(1),'LineWidth',2.5,'LineStyle','-','color',[254/255    157/255    46/255]);
set(hnew(5),'LineWidth',2.5,'LineStyle','-','color',[55/255    226/255    147/255]);
set(hnew(3),'LineWidth',2.5,'LineStyle','-','color',[42/255    83/255    224/255]);
set(hnew(4),'LineWidth',2.5,'LineStyle','-','color',[10/255    10/255    200/255]);

xlabel( 'Age','FontSize',52);
ylabel('Developmental velocity','FontSize',52);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
set(gca,'position',[ 0.25    0.24   0.71  0.71]);
set(gcf, 'PaperPositionMode', 'auto');
legend('off');
print(gcf,'-dtiff','-r300','Global_Slope');
close all


