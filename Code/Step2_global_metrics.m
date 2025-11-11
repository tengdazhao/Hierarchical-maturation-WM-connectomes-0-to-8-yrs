%Global metrics 
% using subset1 as an example
% For multicenter fitting, using multi_center_multimodelnlm_global.m instead of multimodelnlm_global.m
% For multicenter fitting, network atrriobutes should be harmonized (Adj_site_effect.m) before using multi_center_multimodelnlm_global.m
% 

% load graph theoritical measurements
load([Pathnew,'Forgithub\Data\Gretna_results\Subset1\gretna_results_subset1.mat']);
load([Pathnew,'Forgithub\Data\Network\Subset1\Matrixallpar500FN_subset1.mat']);

Matrixallpar500FN=Matrixallpar500FN_subset1;
Age_gender=Age_gender_Brainsize_subset1(:,1:2);
Brainsz=Age_gender_Brainsize_subset1(:,3);
Eg=Eg_subset1;
Eloc=Eloc_subset1;


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

[indminEg mdallEg notEg]=multimodelnlm_global(Age_gender_sp,Brainsz,Eg,'Global Eff',plotornot,indall_select,1)
[indminEloc mdallEloc notEloc]=multimodelnlm_global(Age_gender_sp,Brainsz,Eloc,'Local_E',plotornot,indall_select,1)


%growth rates
Age=Age_gender(indall_select,1)./30./12;
%%%Slope global
h=plot(mdallEg{1,1}{2,indminEg(1,1)},Age,mdallEg{1,1}{1,indminEg(1,1)}.Fitted, 'predfun');
hxvalue(:,1)=h(2,1).XData;hyvalue(:,1)=h(2,1).YData;hyhvalue(:,1)=h(3,1).YData;hylvalue(:,1)=h(4,1).YData;
h=plot(mdallEloc{1,1}{2,indminEloc(1,1)},Age,mdallEloc{1,1}{1,indminEloc(1,1)}.Fitted, 'predfun');hold on;
hxvalue(:,2)=h(2,1).XData;hyvalue(:,2)=h(2,1).YData;hyhvalue(:,2)=h(3,1).YData;hylvalue(:,2)=h(4,1).YData;

close all
Gradeintsall=[];Gradeintshall=[];Gradeintslall=[];
for ii=1:2
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

set(hnewl(2),'LineWidth',0.5,'LineStyle','-','color',[255/255 255/255 255/255]);
set(hnewl(1),'LineWidth',0.5,'LineStyle','-','color',[255/255    255/255    255/255]);
% 
fill([hnewh(2).XData,fliplr(hnewh(2).XData)],[hnewh(2).YData,fliplr(hnewl(2).YData)],[227/255 150/255 120/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色
fill([hnewh(1).XData,fliplr(hnewh(1).XData)],[hnewh(1).YData,fliplr(hnewl(1).YData)],[254/255    210/255    150/255],'EdgeColor','none','FaceAlpha',0.5);%画填充颜色

set(hnew(2),'LineWidth',2.5,'LineStyle','-','color',[227/255 68/255 10/255]);
set(hnew(1),'LineWidth',2.5,'LineStyle','-','color',[254/255    157/255    46/255]);

xlabel( 'Age','FontSize',52);
ylabel('Developmental velocity','FontSize',52);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
set(gca,'position',[ 0.25    0.24   0.71  0.71]);
set(gcf, 'PaperPositionMode', 'auto');
legend('off');
print(gcf,'-dtiff','-r300','Global_Slope');
close all


