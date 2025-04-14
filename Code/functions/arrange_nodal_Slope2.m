function [FourModelsigthre]=arrange_nodal_Slope2(indminNodalStrmdall,indminNodalindre,notenoughNodalStr,Pathnew,opt,scale,Pathnewplot,pthretype,result_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pathnew=Pathnew(1:end-24);
load([Pathnew,'\Data\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);
%load([Pathnew,'\Data\Cov\Age_gender_brainsize.mat']);
md4plot=indminNodalStrmdall{1}{1,1}{1,1};
    tbl=md4plot.Variables;%Using fit to constrain the lower boundary of beta2 above ZERO (Or the model will be failed).
    Age=tbl.Age;
nodal_eff=Ne';
dirsvr=(Pathnewplot);
finalYCorr=mean(nodal_eff(:,find((Age>7)&(Age<=9))),2);

%% re-arranging fitting results
roi_num=length(indminNodalindre);
edges=zeros(roi_num,roi_num);
cd(Pathnewplot);
reompr=[];
for i=1:roi_num
     temres=indminNodalindre{i}.ModelselectRestults;
   modelorder=indminNodalindre{i}.ModelselectRestults(3,:);
    [notus temres(1,:)]=sort(modelorder,'ascend');
    if(temres(1,1)==3)
        if((abs(modelorder(3)-modelorder(2)))<2)
       temres(1,1)=2;temres(1,2)=3;
        end
    end
    
    if(temres(1,1)==2)
        if((abs(modelorder(2)-modelorder(1)))<2)
       temres(1,1)=1;temres(1,2)=2;
        end
    end
    reompr(1,:,i)=temres(1,:);%AICc 
    reompr(2,:,i)=temres(5,temres(1,:));%p
    reompr(3,:,i)=temres(6,temres(1,:));%r
    reompr(4,:,i)=temres(6,temres(1,:))+temres(1,:);%r for plot
    if((reompr(4,:,i)>4)&(reompr(4,:,i)<5))
        reompr(4,:,i)=reompr(4,:,i)+2;
    end
end
reomprtemp=reompr;
reompr4=[];
for i=1:roi_num
reompr4(:,:,i)=reompr(:,find(reomprtemp(1,:,i)<999),i);
end

allp=squeeze([reompr4(2,1,:)]);
allm=squeeze(reompr4(1,1,:));
fdrp=[];[fdrp(1),fdrp(2)]=gretna_FDR(allp,0.01);

if (pthretype==1)
pthref=fdrp(1);%Final FDRp
Pnametreh=['FDR'];
elseif(pthretype==0)
pthref=0.01;%Final uncoorectp
Pnametreh=['Uncorrp001'];
end
nm=[1 2 3];
result_namethre=[result_name,'_',Pnametreh];
reompr4thre1=reompr4;
for i = 1:length(nm)
    indt=find(reompr4(1,1,:)==nm(i));
    nummper(i)=length(indt);
    reompr4thre1(2,1,indt(find(reompr4thre1(2,1,indt)>(pthref))))=999;%0.01
end
% 
indnodes=find(reompr4thre1(2,1,:)<1);
length(indnodes)



intersectionAge=[];Gradeintsall=[];Index2all=[];Mean_rapid_slope=[];;intersectionpointall=[];Median_rapid_slope=[];
for i=1:length(indnodes);
    md4plot=indminNodalStrmdall{indnodes(i)}{1,1}{1,reompr4thre1(1,1,indnodes(i))};
    fitresult=indminNodalStrmdall{indnodes(i)}{1,1}{2,reompr4thre1(1,1,indnodes(i))};
    tbl=md4plot.Variables;%Using fit to constrain the lower boundary of beta2 above ZERO (Or the model will be failed).
    Data_rremain=tbl.data_r-md4plot.Coefficients.Estimate(1)-(md4plot.Coefficients.Estimate(end-1).*tbl.Gender)-(md4plot.Coefficients.Estimate(end).*tbl.BrainSize);
    Data_rremain = Data_rremain+(mean(tbl.data_r)-mean(Data_rremain));
    Age=tbl.Age;
    name=[result_name,'-',num2str(indnodes(i)),'-md',num2str(reompr4thre1(1,1,indnodes(i))),'.tif'];
    h = figure('visible','off');
 %   h = figure;
    h=plot(fitresult,Age,Data_rremain, 'predobs');
    Age2=get(h(2),'xdata');
    Index2=get(h(2),'ydata');
    Dx=diff(Age2);
    Dy=diff(Index2); 
    Dxy=Dy./Dx;
    Gradeintsall(i,:)=Dxy;
    Index2all(i,:)=Index2;
    %figure;plot(fitresult,Age,Data_rremain, 'predobs');
    %hold on;
    %plot(Age2(2:end),Dxy)
    Ypoint=zeros(1,length(Index2))+((8/9).*(max(Index2)-min(Index2)))+min(Index2);
    Fp=(Ypoint)-Index2;
    Fp(find(Fp>0))=1;Fp(find(Fp<0))=-1;
    dffFp=diff(Fp);
    intersectionpointall(i)=find(dffFp~=0);
    intersectionAge(i)=Age2(intersectionpointall(i));
    Mean_rapid_slope(i)=mean(Gradeintsall(i,1:intersectionpointall(i)));
    Median_rapid_slope(i)=mean(Gradeintsall(i,1:intersectionpointall(i)));
    Index2;
    %hold on;
    %plot(intersectionAge(i),Ypoint(1),'*');
end


%%Micro

allmicromask=dir([Pathnew,'\Data\Micro\*.mat']);
ADavg_streamed_node=[];FAavg_streamed_node=[];MDavg_streamed_node=[];RDavg_streamed_node=[];
for i =1:size(allmicromask,1);
load([allmicromask(i).folder,'\',allmicromask(i).name]);
FAavg_streamed_node(i,:) =FAavg_mat_node;
RDavg_streamed_node(i,:) =RDavg_mat_node;
end

mean_FAavg_streamed_node=mean(FAavg_streamed_node,1);
mean_RDavg_streamed_node=mean(RDavg_streamed_node,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sig


result_namethre='Slope_FinalMatu';
load([Pathnew,'\Data\Atlas\Her3name.mat']);
finalYCorrSys=[];GradeintsallsigSys=[];
for i=1:max(AALVonE(:,2));
indsys=find(AALVonE(:,2)==i);
indsyssig=intersect(indnodes,indsys);
finalYCorrSys{1,i}=finalYCorr(indsyssig);
[rowindx, colindy] = find(ismember(indnodes, indsyssig));
GradeintsallsigSys{1,i}=Gradeintsall(rowindx,:);
end


finalYCorrSysAsso=[finalYCorrSys{2}; finalYCorrSys{3}; finalYCorrSys{7};];
finalYCorrSysPri=[finalYCorrSys{1}; finalYCorrSys{4}; finalYCorrSys{5};finalYCorrSys{6}];
GradeintsallsigSysAsso=[GradeintsallsigSys{1,2}; GradeintsallsigSys{1,3}; GradeintsallsigSys{1,7}];
GradeintsallsigSysPri=[GradeintsallsigSys{1,1}; GradeintsallsigSys{1,4}; GradeintsallsigSys{1,5};GradeintsallsigSys{1,6}];
length(GradeintsallsigSysAsso)

rtesta=[];ptesta=[];
for i=1:length(GradeintsallsigSysAsso)
%Gradeintsall2=Gradeintsall(find(typeofstatsig==2),i);
XX=GradeintsallsigSysAsso(:,i);
YY=finalYCorrSysAsso;
[rtesta(i) ptesta(i)]=corr(XX,YY);
end
rtestp=[];ptestp=[];
for i=1:length(GradeintsallsigSysPri)
%Gradeintsall2=Gradeintsall(find(typeofstatsig==2),i);
XX=GradeintsallsigSysPri(:,i);
YY=finalYCorrSysPri;
[rtestp(i) ptestp(i)]=corr(XX,YY);
end

Fitted_Subgroup=[ round(length(Age2)/1100)....
round(length(Age2)/8.9) ....
round(length(Age2)/3.95)  round(length(Age2)/2.63) ....
round(length(Age2)/1.98) round(length(Age2)/1.60)....
round(length(Age2)/1.34) round(length(Age2)/1.01)];
Age2(Fitted_Subgroup)

[maxrass maxrassind]=max(rtesta);


[p1 p2]=gretna_FDR([ptesta ptestp],0.05);

figure
hap=plot([rtesta; rtestp]')
hold on;
ha3=plot([rtesta(find(ptesta>p1))]');
set(hap(1),'LineWidth',4.6,'LineStyle','-','color',[227/255 211/255 132/255]);
set(hap(2),'LineWidth',4.6,'LineStyle','-','color',[168/255     168/255    168/255]);
set(ha3,'LineWidth',4.6,'LineStyle','-','color',[168/255     168/255    168/255]);
xlabel( 'Age','FontSize',45);
ylabel('Correlation Coefficiency','FontSize',45);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 50 920 800]);
set(gca,'position',[ 0.19    0.21   0.72  0.70]);
set(gca,'XTick',Fitted_Subgroup,'XTickLabel',[0 1.0 2.0 3.0 4.0 5.0 6.0 8.0],'FontSize',40);
set(gca,'YTick',[-0.2 0 0.2 0.4 0.6 0.8],'YTickLabel',[-0.2 0 0.2 0.4 0.6 0.8],'FontSize',40);
set(gca,'fontsize',40,'LineWidth',5,'FontWeight','bold');
set(gcf, 'PaperPositionMode', 'auto');
plot(maxrassind,maxrass,'.','Color','r','MarkerSize',45);
box off
Age2(maxrassind)
print(gcf,'-dtiff','-r300',[result_namethre 'Corr.tif']); close all;



%%% using mean fa and max slope for mediation
[maxrass maxrassind]=max(rtesta);
statsallmed=[];rcall=[];
ii=length(maxrassind)
[Agenew Ageindsort]=sort(Age,'descend');
[nnn slideindex]=min(abs(Agenew-Age2(maxrassind(ii))));

if(slideindex>107)
    Slideindexraw=Ageindsort([slideindex-7:end]);
elseif(slideindex<7)
    Slideindexraw=Ageindsort([1:slideindex+7]);
else
Slideindexraw=Ageindsort([slideindex-7:slideindex+7]);
end
%Mean_forMedation_FA=mean(FAavg_masked_node(Slideindexraw,:));

Mean_forMedation_FA=[];
Mean_forMedation_FA=mean(FAavg_streamed_node(Slideindexraw,:));
MeanAll_forMedation_FASys=[];

FArcall=[];RDrcall=[];
statsallmedFA=[];statsallmedRD=[];
for i=1:max(AALVonE(:,2));
indsys=find(AALVonE(:,2)==i);
indsyssig=intersect(indnodes,indsys);
%MeanAll_forMedation_FASys{1,i}=mean_FAavg_streamed_node(indsyssig)';
MeanAll_forMedation_FASys{1,i}=Mean_forMedation_FA(indsyssig)';

end

MeanAll_forMedation_FASys=[MeanAll_forMedation_FASys{2}; MeanAll_forMedation_FASys{3}; MeanAll_forMedation_FASys{7};];

Age2(maxrassind)
%XX=GradeintsallsigSysAsso(:,maxrassind(ii));
XX=GradeintsallsigSysAsso(:,maxrassind);
YY=finalYCorrSysAsso;
MM=MeanAll_forMedation_FASys;

[paths,statsallmedFA{1,ii}] = mediation(XX,YY,MM,'plot','boot','bootsamples',10000,'doCIs','verbose');
[r1 p1]=corr(XX,MM);
[r2 p2]=partialcorr(MM,YY,XX);
[r3 p3]=corr(XX,YY);
[r4 p4]=partialcorr(XX,YY,MM);
[r5 p5]=corr(YY,MM);
ill=1;
FArcall(1,ill:ill+1,ii)=[r1 p1];
FArcall(2,ill:ill+1,ii)=[r2 p2];
FArcall(3,ill:ill+1,ii)=[r3 p3];
FArcall(4,ill:ill+1,ii)=[r4 p4];
FArcall(5,ill:ill+1,ii)=[statsallmedFA{1,ii}.ste(1,5) statsallmedFA{1,ii}.p(1,5)];
FArcall(6,ill:ill+1,ii)=[statsallmedFA{1,ii}.ci(1,5,1) statsallmedFA{1,ii}.ci(1,5,2)];


[corr_stats]=corrstat(XX,YY,[]);
corrplot(XX,corr_stats,'Nodal growth rates at 3.7 y','Nodal efficiency at 8 y',[],[],[],[],[],['Median_corr_FA.tif'])

[corr_stats]=corrstat(XX,MM,[]);
corrplot(XX,corr_stats,'Nodal growth rates at 3.7 y','Nodal FA at 3.7 y',[],[],[],[],[],['Median_corr2_FA.tif'])


[corr_stats]=corrstat(MM,YY,[]);
corrplot(MM,corr_stats,'Nodal FA at 3.7 y','Nodal efficiency at 8 y',[],[],[],[],[],['Median_corr3_FA.tif'])


Mean_forMedation_RD=[];
Mean_forMedation_RD=mean(RDavg_streamed_node(Slideindexraw,:));
MeanAll_forMedation_RDSys=[];
for i=1:max(AALVonE(:,2));
indsys=find(AALVonE(:,2)==i);
indsyssig=intersect(indnodes,indsys);
%MeanAll_forMedation_FASys{1,i}=mean_FAavg_streamed_node(indsyssig)';
MeanAll_forMedation_RDSys{1,i}=Mean_forMedation_RD(indsyssig)';

end

MeanAll_forMedation_RDSys=[MeanAll_forMedation_RDSys{2}; MeanAll_forMedation_RDSys{3}; MeanAll_forMedation_RDSys{7};];

Age2(maxrassind)
%XX=GradeintsallsigSysAsso(:,maxrassind(ii));
XX=GradeintsallsigSysAsso(:,maxrassind);
YY=finalYCorrSysAsso;
MM=MeanAll_forMedation_RDSys;
[paths,statsallmedRD{1,ii}] = mediation(XX,YY,MM,'plot','boot','bootsamples',10000,'doCIs','verbose');
[r1 p1]=corr(XX,MM);
[r2 p2]=partialcorr(MM,YY,XX);
[r3 p3]=corr(XX,YY);
[r4 p4]=partialcorr(XX,YY,MM);
[r5 p5]=corr(YY,MM);
ill=1;
RDrcall(1,ill:ill+1,ii)=[r1 p1];
RDrcall(2,ill:ill+1,ii)=[r2 p2];
RDrcall(3,ill:ill+1,ii)=[r3 p3];
RDrcall(4,ill:ill+1,ii)=[r4 p4];
RDrcall(5,ill:ill+1,ii)=[statsallmedRD{1,ii}.ste(1,5) statsallmedRD{1,ii}.p(1,5)];
RDrcall(6,ill:ill+1,ii)=[statsallmedRD{1,ii}.ci(1,5,1) statsallmedRD{1,ii}.ci(1,5,2)];


[corr_stats]=corrstat(XX,YY,[]);
corrplot(XX,corr_stats,'Nodal growth rates at 3.7 y','Nodal efficiency at 8 y',[],[],[],[],[],['Median_corr_RD.tif'])

[corr_stats]=corrstat(XX,MM,[]);
corrplot(XX,corr_stats,'Nodal growth rates at 3.7 y','Nodal RD at 3.7 y',[],[],[],[],[],['Median_corr2_RD.tif'])


[corr_stats]=corrstat(MM,YY,[]);
corrplot(MM,corr_stats,'Nodal RD at 3.7 y','Nodal efficiency at 8 y',[],[],[],[],[],['Median_corr3_RD.tif'])



