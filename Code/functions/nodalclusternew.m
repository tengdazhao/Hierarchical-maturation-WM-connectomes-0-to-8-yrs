function [maxpar maxind silhAll idxAll]=nodalclusternew(nodalEg_All_Node,Distmetric,Startmethod,numclu,Reptimes,fname);
Allpar=[];silhAll=[];idxAll=[];
for i=2:numclu
    idx=[];
    for j=1:1
rmpath('H:\Dpan\Program Files\MATLAB\R2011b\toolbox\stats\');
rmpath('H:\Dpan\Program Files\MATLAB\R2011b\toolbox\stats\stats');
rmpath('D:\Program Files\MATLAB\R2022b\NeuroImage_Toolbox\SPM12\external\fieldtrip\external\stats');
rmpath('D:\Program Files\MATLAB\R2020b\toolbox\NeuroImage_Toolbox\CanlabCore-master');
idx(j,:)= kmeans(nodalEg_All_Node',i, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
Allpar(:,i-1)=idx(j,:) ;
figure
[silh(j,:),h] = silhouette(nodalEg_All_Node',idx(j,:),Distmetric);
tryn=1;
silhtryn=[];
while(min(silh(j,:))<-0.05)
[silhtryn(tryn,:),h] = silhouette(nodalEg_All_Node',idx(j,:),Distmetric);
    silh(j,:)=silhtryn(tryn,:);
tryn=tryn+1
if(tryn>10)
    [myv mintry]=max(mean(silhtryn'));
    silh(j,:)=silhtryn(mintry,:);
    break
end
end
    end
    [cn cn2 ]=min(mean(silh'));
idxAll(:,i-1)=idx(cn2,:);
silhAll(:,i-1)=silh(cn2,:);
silhAll(:,i-1)=silh(cn2,:);
set(get(gca,'Children'),'FaceColor',[.8 .8 1])
xlabel('Silhouette Value')
ylabel('Cluster')
print(gcf,'-dtiff','-r500',['sqEuclideanuniformbor.tif']);  
close all
end


% 
% idx3 = kmeans(nodalEg_All_Node',3, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,2)=idx3;
% figure
% [silh3,h] = silhouette(nodalEg_All_Node',idx3,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean3sample.tif']);  
% 
% 
% 
% idx4 = kmeans(nodalEg_All_Node',4,'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,3)=idx4;
% figure
% [silh4,h] = silhouette(nodalEg_All_Node',idx4,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean4.tif']);  
% 
% 
% idx5 = kmeans(nodalEg_All_Node',5, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,4)=idx5;
% figure
% [silh5,h] = silhouette(nodalEg_All_Node',idx5,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean5cluster.tif']);  
% 
% 
% idx6 = kmeans(nodalEg_All_Node',6, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,5)=idx6;
% figure
% [silh6,h] = silhouette(nodalEg_All_Node',idx6,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean6.tif']);  
% 
% idx7 = kmeans(nodalEg_All_Node',7, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,6)=idx7;
% figure
% [silh7,h] = silhouette(nodalEg_All_Node',idx7,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean7.tif']);  
% 
% idx8 = kmeans(nodalEg_All_Node',8, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,7)=idx8;
% figure
% [silh8,h] = silhouette(nodalEg_All_Node',idx8,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean8.tif']);  
% 
% idx9 = kmeans(nodalEg_All_Node',7, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,8)=idx9;
% figure
% [silh9,h] = silhouette(nodalEg_All_Node',idx9,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean9.tif']);  
% 
% idx10 = kmeans(nodalEg_All_Node',10, 'dist',Distmetric, 'Start',Startmethod,'display','iter','OnlinePhase','on','Replicates',Reptimes);
% Allpar(:,9)=idx10;
% figure
% [silh10,h] = silhouette(nodalEg_All_Node',idx10,Distmetric);
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])
% xlabel('Silhouette Value')
% ylabel('Cluster')
% print(gcf,'-dtiff','-r500',['sqEuclidean10.tif']);  


barslope = [mean(silhAll); std(silhAll);];
[nn maxind]=max(barslope(1,:));
maxpar=Allpar(:,maxind);
figure;
%sca = [0.5 0.2];
axis([0 length(barslope)+1 min(mean(silhAll))-0.05 max(mean(silhAll))+0.2]);
hold on;
set(gcf,'position', [69          69        1091         992]);
set(gca,'position',[0.18,0.18,0.75,0.75]);
set(gcf, 'PaperPositionMode', 'auto')  
set(gca, 'fontsize', 38,'LineWidth',4,'FontWeight','bold');

% set(gca,'YTick',[0.1 0.3 0.5]);
% x_tick = {'TT','GG+GT'}
% set(gca, 'XtickLabel',x_tick);
ylabel('Silhouette Values');
mvb = barslope(1,:)
%mvb(2,:)= NaN;
bara = bar(mvb','grouped','LineWidth',0.4,'facecolor',[0.9 0 0],'BarWidth',0.68);
set(bara,'FaceColor',[212/255 208/255 200/255]);
% set(bara(2),'FaceColor',[212/255 208/255 200/255]);
% set(bara(3),'FaceColor',[212/255 208/255 200/255]);
% set(bara(4),'FaceColor',[212/255 208/255 200/255]);
% set(bara(5),'FaceColor',[212/255 208/255 200/255]);

x= [1:1:numclu-1];
%errorbar_ztd(barslope(1,:)',barslope(2,:)',x,0.32,4);
set(gca,'XTick',x,'XTicklabel',[' 2'; ' 3'; ' 4'; ' 5';]);
set(gcf,'Color',[1 1 1]) ;
%legend('2 Clusters','3 Clusters','4 Clusters');
% set(legend,'Position',[ 0.66    0.74   0.32    0.23],'fontsize', 30,'FontWeight','bold');
% text(0.63,max(sca)-(max(sca)-min(sca))/120, 'p = 6.5 x 10^{\fontsize{28}-4}','fontsize', 40,'FontWeight','bold');
print(gcf,'-dtiff','-r300',fname);
