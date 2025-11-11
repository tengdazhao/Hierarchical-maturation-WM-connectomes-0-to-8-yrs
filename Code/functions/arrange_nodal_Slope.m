function [FourModelsigthre]=arrange_nodal_Slope(Nodal_metric,indminNodalStrmdall,indminNodalindre,notenoughNodalStr,Pathnew,opt,scale,Pathnewplot,pthretype,result_name)
%reorder
roi_num=length(indminNodalindre);
edges=zeros(roi_num,roi_num);
cd(Pathnewplot);
reompr=[];
for i=1:roi_num
     temres=indminNodalindre{i}.ModelselectRestults;
    %temres(1,:)=[notenoughNodalStr{i}.select setdiff(temres(1,:),notenoughNodalStr{i}.select)];
   modelorder=indminNodalindre{i}.ModelselectRestults(3,:);
   %temres(1,:)=modelorder;
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

%not use
reomprtemp=reompr;
% reomprtemp(1,find(reompr(1,:,:)==5))=999;%#ok<FNDSB>
% reomprtemp(1,find(reompr(1,:,:)==6))=999;%#ok<FNDSB>
% reomprtemp(1,find(reompr(1,:,:)==2))=999;%#ok<FNDSB>
% 
% 
reompr4=[];
for i=1:roi_num
reompr4(:,:,i)=reompr(:,find(reomprtemp(1,:,i)<999),i);
end
% 
% for ii=1:size(reompr4,3)
% reompr4(2,find(reompr4(3,:,ii)<0.11),ii)=999;%R>0.05
% end

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
% FourModelsigthre{1,1}=length(indnodes);
% FourModelsigthre{1,2}=indnodes;
% FourModelsigthre{1,3}=reompr4thre1(1,:,indnodes);
load([Pathnew(1:end-14),'\Her3name.mat']);
%load([Pathnew(1:end-14),'\Her3nameR2.mat']);
PerHer4=[];
for i=1:max(AALHer4(:,2))
    tempm=allm(find(AALHer4(:,2)==i));   
    PerHer4(i,1)=length(tempm(find(tempm==1)))./length(tempm);
    PerHer4(i,2)=length(tempm(find(tempm==2)))./length(tempm);
end
PerVonE=[];
for i=1:max(AALVonE(:,2))
    tempm=allm(find(AALVonE(:,2)==i));
    PerVonE(i,1)=length(tempm(find(tempm==1)))./length(tempm);
    PerVonE(i,2)=length(tempm(find(tempm==2)))./length(tempm);
end
PerYeo=[];
for i=1:max(AALYeo(:,2))
    tempm=allm(find(AALYeo(:,2)==i));
    PerYeo(i,1)=length(tempm(find(tempm==1)))./length(tempm);
    PerYeo(i,2)=length(tempm(find(tempm==2)))./length(tempm);
end
% 
% bar(PerHer4(:,1)+PerHer4(:,2),'FaceColor',[237/255 176/255 33/255]);hold on;bar(PerHer4(:,1),'FaceColor',[70/255 171/255,215/255]);
% bar(PerVonE(:,1)+PerVonE(:,2),'FaceColor',[237/255 176/255 33/255]);hold on;bar(PerVonE(:,1),'FaceColor',[70/255 171/255,215/255]);
% bar(PerYeo(:,1)+PerYeo(:,2),'FaceColor',[237/255 176/255 33/255]);hold on;bar(PerYeo(:,1),'FaceColor',[70/255 171/255,215/255]);
% 
%pthref=0.05;%Final FDRp
PerHer4=[]; Forchi=[];
for i=1:max(AALHer4(:,2))
    tempm=allm(find(AALHer4(:,2)==i));  
    tempp=allp(find(AALHer4(:,2)==i));  
    tempmSig=tempm(find(tempp<=pthref)); 
    PerHer4(i,1)=length(tempmSig(find(tempmSig==1)))./length(tempm);
    PerHer4(i,2)=length(tempmSig(find(tempmSig==2)))./length(tempm);
    PerHer4(i,3)=(length(tempm)-length(tempmSig))./length(tempm);
    Forchi{i}=[length(tempmSig(find(tempmSig==1))) length(tempmSig(find(tempmSig==2))) ];
end

pHer4=[];QHer4=[];
for  it=1:4;
    for jt=1:4;
        if(it~=jt)
    [pHer4(it,jt), QHer4(it,jt)]= chi2test([Forchi{1,it};Forchi{1,jt}]);
        else
          pHer4(it,jt)=999;  QHer4(it,jt)=999
        end
    end
end



PerVonE=[];ForchiVonE=[];
for i=1:max(AALVonE(:,2))
    tempm=allm(find(AALVonE(:,2)==i));  
    tempp=allp(find(AALVonE(:,2)==i));  
    tempmSig=tempm(find(tempp<=pthref)); 
    PerVonE(i,1)=length(tempmSig(find(tempmSig==1)))./length(tempm);
    PerVonE(i,2)=length(tempmSig(find(tempmSig==2)))./length(tempm);
    PerVonE(i,3)=(length(tempm)-length(tempmSig))./length(tempm);
    ForchiVonE{i}=[length(tempmSig(find(tempmSig==1))) length(tempmSig(find(tempmSig==2))) ];
end


pVonE=[];QVonE=[];
for  it=1:7;
    for jt=1:7;
        if(it~=jt)
    [pVonE(it,jt), QVonE(it,jt)]= chi2test([ForchiVonE{1,it};ForchiVonE{1,jt}]);
        else
          pVonE(it,jt)=999;  QVonE(it,jt)=999
        end
    end
end


PerYeo=[];ForchiYeo=[];
for i=1:max(AALYeo(:,2))
    tempm=allm(find(AALYeo(:,2)==i));  
    tempp=allp(find(AALYeo(:,2)==i));  
    tempmSig=tempm(find(tempp<=pthref)); 
    PerYeo(i,1)=length(tempmSig(find(tempmSig==1)))./length(tempm);
    PerYeo(i,2)=length(tempmSig(find(tempmSig==2)))./length(tempm);
    PerYeo(i,3)=(length(tempm)-length(tempmSig))./length(tempm);
    ForchiYeo{i}=[length(tempmSig(find(tempmSig==1))) length(tempmSig(find(tempmSig==2))) ];
end

pVonYeo=[];QVonYeo=[];
for  it=1:7;
    for jt=1:7;
        if(it~=jt)
    [pVonYeo(it,jt), QVonYeo(it,jt)]= chi2test([ForchiYeo{1,it};ForchiYeo{1,jt}]);
        else
          pVonYeo(it,jt)=999;  QVonYeo(it,jt)=999
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 classes
figure;
set(gcf,'position', [4         357       524         350]);
set(gca,'position',[0.42,0.18,0.55,0.65]);
[xx numone]=sort(PerHer4(:,1),'ascend');
PerHer4new=PerHer4(numone,:);
%barh(PerHer4new(:,1)+PerHer4new(:,2)+PerHer4new(:,3),'FaceColor',[70/255 171/255,215/255]);hold on;barh(PerHer4new(:,2)+PerHer4new(:,3),'FaceColor',[237/255 176/255 33/255]);hold on;barh(PerHer4new(:,3),'FaceColor',[152/255 152/255 152/255]);
barh(PerHer4new(:,1)+PerHer4new(:,2)+PerHer4new(:,3),'FaceColor',[152/255 152/255 152/255]);hold on;barh(PerHer4new(:,1)+PerHer4new(:,2),'FaceColor',[237/255 176/255 33/255]);hold on;barh(PerHer4new(:,1),'FaceColor',[70/255 171/255,215/255]);
set(gca,'YTick',(1:1:4),'YTickLabel',[sysname(1,numone)]);
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
print(gcf,'-dtiff','-r300',[result_namethre 'PerHer4Num.tif']); close all;

[xx numone]=sort(PerHer4(:,1),'descend');
%Her4numonefinal=flip(numone);
pHer4new=pHer4(numone,numone);
figure;
pHer4plot=-log(pHer4new);
pHer4plot(find(pHer4new==999))=0;
pall=[];
pall=pHer4new(intersect((find(triu(pHer4new)>0)),(find(triu(pHer4new)<1))));
[pmmfdr1 pmmfdr2]=gretna_FDR(pall,0.05);
pmmfdr2=0.05;
pHer4plot(find(pHer4new>(pmmfdr2)))=0;
% 
data=(pHer4plot);
% 创建一个绘图函数，只在圆内填充伪彩色的颜色
figure;
hold on;
[row, col] = size(data);
[x, y] = meshgrid(1:col, 1:row);
radius = 0.4; % 圆的半径

% 创建一个颜色映射
cmap = nclCM(365);
cmap(1,:)=[0.89 0.89 0.89];
for i = 1:row
    for j = 1:col
        % 计算圆心的坐标
        cx = j;
        cy = i;
        
        % 计算颜色索引，限制在1到256之间
        color_idx = round(data(i, j)./max(max(data)) * (length(cmap)-1))+1;
        
        % 从颜色映射中获取颜色
        color = cmap(color_idx, :);
        
        rectangle('Position', [j - radius, i - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', color, 'EdgeColor', 'none');
    end
end
axis equal;
colormap(cmap);
colorbar;
axis ij %// use if you want the origin at the top, like with imagesc
axis off
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
hold off;
print(gcf,'-dtiff','-r300',[result_namethre 'NoPerHer4PCirclemap.tif']); close all;

% 7 classes strucural
figure;
set(gcf,'position', [4         357       524         350]);
set(gca,'position',[0.42,0.18,0.55,0.65]);
[xx numone]=sort(PerVonE(:,1),'ascend');
PerVonEnew=PerVonE(numone,:);
llw=1;
barh(PerVonEnew(:,1)+PerVonEnew(:,2)+PerVonEnew(:,3),'FaceColor',[0.89 0.89 0.89],'LineWidth',llw);hold on;barh(PerVonEnew(:,1)+PerVonEnew(:,2),'FaceColor',[237/255 176/255 33/255],'LineWidth',llw);hold on;barh(PerVonEnew(:,1),'FaceColor',[70/255 171/255,215/255],'LineWidth',llw);
set(gca,'YTick',(1:1:7),'YTickLabel',[sysname(2,numone)]);
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'FontWeight','bold');
%axis off
print(gcf,'-dtiff','-r300',[result_namethre 'PerVonENumtif']); close all;

[xx numone]=sort(PerVonE(:,1),'descend');
%VonEnumonefinal=flip(numone);
pVonEnew=pVonE(numone,numone);
pVonEplot=-log(pVonEnew);
pVonEplot(find(pVonEnew==999))=0;
pall=[];
pall=pVonEnew(intersect((find(triu(pVonEnew)>0)),(find(triu(pVonEnew)<1))));
[pmmfdr1 pmmfdr2]=gretna_FDR(pall,0.05);
pVonEplot(find(pVonEnew>(pmmfdr1)))=0;


data=(pVonEplot);
% 创建一个绘图函数，只在圆内填充伪彩色的颜色
figure;
hold on;
[row, col] = size(data);
[x, y] = meshgrid(1:col, 1:row);
radius = 0.4; % 圆的半径

% 创建一个颜色映射
cmap = nclCM(365);
cmap(1,:)=[0.89 0.89 0.89];
for i = 1:row
    for j = 1:col
        % 计算圆心的坐标
        cx = j;
        cy = i;
        
        % 计算颜色索引，限制在1到256之间
        color_idx = round(data(i, j)./max(max(data)) * (length(cmap)-1))+1;
        
        % 从颜色映射中获取颜色
        color = cmap(color_idx, :);
        
        rectangle('Position', [j - radius, i - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', color, 'EdgeColor', 'none');
    end
end

axis equal;
colormap(cmap);
colorbar;
axis ij %// use if you want the origin at the top, like with imagesc
axis off
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
hold off;
print(gcf,'-dtiff','-r300',[result_namethre 'PerVonEPCirclemap.tif']); close all;

% 7 classes functional
figure;
set(gcf,'position', [4         357       524         350]);
set(gca,'position',[0.42,0.18,0.55,0.65]);
[xx numone]=sort(PerYeo(:,1),'ascend');
PerYeonew=PerYeo(numone,:);
barh(PerYeonew(:,1)+PerYeonew(:,2)+PerYeonew(:,3),'FaceColor',[152/255 152/255 152/255]);hold on;barh(PerYeonew(:,1)+PerYeonew(:,2),'FaceColor',[237/255 176/255 33/255]);hold on;barh(PerYeonew(:,1),'FaceColor',[70/255 171/255,215/255]);
set(gca,'YTick',(1:1:7),'YTickLabel',[sysname(3,numone)])
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
print(gcf,'-dtiff','-r300',[result_namethre '_PerYeoNum.tif']); close all;

%Yeonumonefinal=flip(numone);
[xx numone]=sort(PerYeo(:,1),'descend');
pVonYeonew=pVonYeo(numone,numone);
figure;
pYeoEplot=-log(pVonYeonew);
pYeoEplot(find(pVonYeonew==999))=0;
pall=[];
pall=pVonYeonew(intersect((find(triu(pVonYeonew)>0)),(find(triu(pVonYeonew)<1))));
[pmmfdr1 pmmfdr2]=gretna_FDR(pall,0.05);
pmmfdr2=0.05;
pYeoEplot(find(pVonYeonew>(pmmfdr2)))=0;
% 

data=(pYeoEplot);
% 创建一个绘图函数，只在圆内填充伪彩色的颜色
figure;
hold on;
[row, col] = size(data);
[x, y] = meshgrid(1:col, 1:row);
radius = 0.4; % 圆的半径

% 创建一个颜色映射
cmap = nclCM(365);
cmap(1,:)=[0.89 0.89 0.89];
for i = 1:row
    for j = 1:col
        % 计算圆心的坐标
        cx = j;
        cy = i;
        
        % 计算颜色索引，限制在1到256之间
        color_idx = round(data(i, j)./max(max(data)) * (length(cmap)-1))+1;
        
        % 从颜色映射中获取颜色
        color = cmap(color_idx, :);
        
        rectangle('Position', [j - radius, i - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', color, 'EdgeColor', 'none');
    end
end

axis equal;
colormap(cmap);
colorbar;
axis ij %// use if you want the origin at the top, like with imagesc
axis off
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
hold off;
print(gcf,'-dtiff','-r300',[result_namethre 'NoPerYeoPCirclemap.tif']); close all;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%slope stat

intersectionAge=[];Gradeintsall=[];Index2all=[];Mean_rapid_slope=[];;intersectionpointall=[];Median_rapid_slope=[];
for i=1:length(indnodes);
    try
    md4plot=indminNodalStrmdall{indnodes(i)}{1,1}{1,reompr4thre1(1,1,indnodes(i))};
    fitresult=indminNodalStrmdall{indnodes(i)}{1,1}{2,reompr4thre1(1,1,indnodes(i))};
    catch
    md4plot=indminNodalStrmdall{indnodes(i)}{1,reompr4thre1(1,1,indnodes(i))};
    fitresult=indminNodalStrmdall{indnodes(i)}{2,reompr4thre1(1,1,indnodes(i))};
    end
    
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
    difftemp=find(dffFp~=0);
    intersectionpointall(i)=difftemp(1);
    intersectionAge(i)=Age2(intersectionpointall(i));
    Mean_rapid_slope(i)=mean(Gradeintsall(i,1:intersectionpointall(i)));
    Median_rapid_slope(i)=mean(Gradeintsall(i,1:intersectionpointall(i)));
    Index2;
    %hold on;
    %plot(intersectionAge(i),Ypoint(1),'*');
end

typeofstat=zeros(roi_num,1);
typeofstat(indnodes)=reompr4(1,1,indnodes);

Fitted_Subgroup=[ round(length(Age2)/22)  round(length(Age2)/8.9) round(length(Age2)/3.95) round(length(Age2)/1.98) round(length(Age2)/1.15)];
Age2(Fitted_Subgroup)

Fitted_Subgroup=[ round(length(Age2)/1100)....
    round(length(Age2)/21)  round(length(Age2)/8.9) round(length(Age2)/5.6)....
    round(length(Age2)/3.95) round(length(Age2)/3.15)  round(length(Age2)/2.63) round(length(Age2)/2.23)....
    round(length(Age2)/1.98) round(length(Age2)/1.60)....
    round(length(Age2)/1.34) round(length(Age2)/1.01)];
Age2(Fitted_Subgroup)

%Fitted_Subgroup=[ round(length(Age2)/1100) round(length(Age2)/3.95) round(length(Age2)/2.63) ....
%    round(length(Age2)/1.98) round(length(Age2)/1.60)....
%    round(length(Age2)/1.01)];

sizenode=zeros(roi_num,1);
colornode=zeros(roi_num,1);
tempcolornode=zeros(roi_num,1);
Gradeintsallplot=Gradeintsall;
%opt=[Pathnewplot,'\Slpoptnew.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%slope brain nodes stats
load(['Her3name.mat']);
InterAgeSys=[];
for i=1:max(AALVonE(:,2))
    indsyssiginnodes=[];
    indsys=find(AALVonE(:,2)==i);  
    indsyssig=intersect(indnodes,indsys);
    for iii=1:length(indsyssig);indsyssiginnodes(iii)=find(indnodes==indsyssig(iii));end
    InterAgeSys{1,i}=intersectionAge(indsyssiginnodes)';
end


InterAgeSysAsso=[InterAgeSys{2}; InterAgeSys{3}; InterAgeSys{7};];
InterAgeSysPri=[InterAgeSys{1}; InterAgeSys{4}; InterAgeSys{5};InterAgeSys{6}];

%slope brain nodes plot

sizenode(indnodes)=reompr4thre1(3,1,indnodes);
colornode(indnodes)=intersectionAge;
%plotnode(colornode,sizenode,edges,Pathnew,opt,scale,Pathnewplot,[result_namethre,'MaturationAge']);
sizenode(indnodes)=reompr4thre1(3,1,indnodes);
colornode(indnodes)=Mean_rapid_slope;
plotnode(colornode,sizenode,edges,Pathnew,opt,scale,Pathnewplot,[result_namethre,'mean_slope']);
sizenode(indnodes)=reompr4thre1(3,1,indnodes);
colornode(indnodes)=Median_rapid_slope;
%plotnode(colornode,sizenode,edges,Pathnew,opt,scale,Pathnewplot,[result_namethre,'mediation_slope']);
for ii=1:length(Fitted_Subgroup)
%for ii=4
colornode(indnodes)=Gradeintsallplot(:,Fitted_Subgroup(ii));
tempcolornode(indnodes,ii)=Gradeintsallplot(:,Fitted_Subgroup(ii));
result_namethre=[Pnametreh,'threFDRSlp',num2str(ii)];
plotnode(tempcolornode(:,ii),sizenode,edges,Pathnew,opt,scale,Pathnewplot,[result_namethre,'_Slope_Agenew',num2str(Age2(Fitted_Subgroup(ii)))]);
close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slopeallmatrix
myslopecmap=[0 0 0.515625;0 0 0.552094256025028;0 0 0.588563512050057;0 0 0.625032768075085;0 0 0.661502024100113;0 0 0.697971280125142;0 0 0.73444053615017;0 0 0.770909792175199;0 0 0.807379048200227;0 0 0.843848304225255;0 0 0.880317560250284;0 0 0.916786816275312;0 0 0.95325607230034;0 0.00628502423763195 0.983440304087737;0 0.0285050022346753 0.997689582115722;0 0.0626638403754253 1;0 0.0991330964004536 1;0 0.135602352425482 1;0 0.17207160845051 1;0 0.208540864475538 1;0 0.245010120500567 1;0 0.281479376525595 1;0 0.317948632550623 1;0 0.354417888575652 1;0 0.39088714460068 1;0 0.427356400625708 1;0 0.463825656650737 1;0 0.500294912675765 1;0 0.536764168700793 1;0 0.577204270808264 1;0 0.62282440285478 1;0 0.649993782419262 1;0 0.661673363719414 1;0 0.673352945019567 1;0 0.685032526319719 1;0 0.696712107619871 1;0 0.708391688920024 1;0 0.720071270220176 1;0 0.731750851520328 1;0 0.743430432820481 1;0 0.755110014120633 1;0 0.766789595420785 1;0 0.778469176720938 1;0 0.79014875802109 1;0 0.801828339321242 1;0 0.813507920621395 1;0 0.825187501921547 1;0 0.836867083221699 1;0 0.848546664521852 1;0 0.860226245822004 1;0 0.871905827122156 1;0 0.883585408422309 1;0 0.895264989722461 1;0 0.906944571022614 1;0 0.918624152322766 1;0 0.930303733622918 1;0 0.941983314923071 1;0 0.953662896223223 1;0 0.965907618554027 1;0 0.976852803387951 1;0 0.987797988221875 1;0.00265843048684292 0.996084742568956 0.997341569513157;0.00990367300120919 0.999784684888513 0.990096326998791;0.0206335427236458 1 0.979366457276354;0.0315787275575695 1 0.96842127244243;0.0425239123914932 1 0.957476087608507;0.0534690972254169 1 0.946530902774583;0.0644142820593406 1 0.935585717940659;0.0753594668932643 1 0.924640533106736;0.0863046517271879 1 0.913695348272812;0.0972498365611115 1 0.902750163438888;0.108195021395035 1 0.891804978604965;0.119140206228959 1 0.880859793771041;0.130085391062882 1 0.869914608937118;0.141030575896806 1 0.858969424103194;0.15197576073073 1 0.848024239269271;0.162920945564653 1 0.837079054435347;0.173866130398577 1 0.826133869601423;0.1848113152325 1 0.8151886847675;0.195756500066424 1 0.804243499933576;0.206701684900348 1 0.793298315099652;0.217646869734271 1 0.782353130265729;0.228592054568195 1 0.771407945431805;0.239537239402118 1 0.760462760597882;0.250482424236042 1 0.749517575763958;0.261427609069966 1 0.738572390930034;0.272372793903889 1 0.727627206096111;0.283317978737813 1 0.716682021262187;0.294263163571736 1 0.705736836428264;0.30520834840566 1 0.69479165159434;0.316153533239584 1 0.683846466760416;0.327098718073507 1 0.672901281926493;0.338043902907431 1 0.661956097092569;0.348989087741354 1 0.651010912258646;0.359934272575278 1 0.640065727424722;0.370879457409202 1 0.629120542590798;0.381824642243125 1 0.618175357756875;0.392769827077049 1 0.607230172922951;0.403715011910972 1 0.596284988089028;0.414660196744896 1 0.585339803255104;0.42560538157882 1 0.57439461842118;0.436550566412743 1 0.563449433587257;0.447495751246667 1 0.552504248753333;0.45844093608059 1 0.54155906391941;0.469386120914514 1 0.530613879085486;0.480331305748438 1 0.519668694251562;0.491276490582361 1 0.508723509417639;0.502221675416285 1 0.497778324583715;0.513166860250209 1 0.486833139749791;0.524112045084132 1 0.475887954915868;0.535057229918056 1 0.464942770081944;0.54600241475198 1 0.45399758524802;0.556947599585904 1 0.443052400414096;0.567892784419828 1 0.432107215580172;0.578837969253751 1 0.421162030746249;0.589783154087675 1 0.410216845912325;0.600728338921599 1 0.399271661078401;0.611673523755523 1 0.388326476244477;0.622618708589447 1 0.377381291410553;0.63356389342337 1 0.36643610657663;0.644509078257294 1 0.355490921742706;0.655454263091218 1 0.344545736908782;0.666399447925142 1 0.333600552074858;0.677344632759066 1 0.322655367240934;0.688289817592989 1 0.311710182407011;0.699235002426913 1 0.300764997573087;0.710180187260837 1 0.289819812739163;0.721125372094761 1 0.278874627905239;0.732070556928685 1 0.267929443071315;0.743015741762608 1 0.256984258237392;0.753960926596532 1 0.246039073403468;0.764906111430456 1 0.235093888569544;0.77585129626438 1 0.22414870373562;0.786796481098303 1 0.213203518901696;0.797741665932227 1 0.202258334067773;0.808686850766151 1 0.191313149233849;0.819632035600075 1 0.180367964399925;0.830577220433999 1 0.169422779566001;0.841522405267922 1 0.158477594732078;0.852467590101846 1 0.147532409898154;0.86341277493577 1 0.13658722506423;0.874357959769694 1 0.125642040230306;0.885303144603618 1 0.114696855396382;0.896248329437542 1 0.103751670562458;0.907193514271465 1 0.0928064857285346;0.918138699105389 1 0.0818613008946108;0.929083883939313 1 0.070916116060687;0.940029068773237 1 0.0599709312267632;0.950974253607161 1 0.0490257463928394;0.961919438441084 1 0.0380805615589156;0.972864623275008 1 0.0271353767249918;0.983809808108932 1 0.016190191891068;0.994237052812132 0.999482059869276 0.0057629471878684;0.999392283579947 0.993692105803168 0.000607716420052767;1 0.983354637389297 0;1 0.972409452555373 0;1 0.961464267721449 0;1 0.950519082887525 0;1 0.939573898053601 0;1 0.928628713219678 0;1 0.917683528385754 0;1 0.90673834355183 0;1 0.895793158717906 0;1 0.884847973883983 0;1 0.873902789050059 0;1 0.862957604216135 0;1 0.852012419382211 0;1 0.841067234548287 0;1 0.830122049714364 0;1 0.81917686488044 0;1 0.808231680046516 0;1 0.797286495212592 0;1 0.786341310378668 0;1 0.775396125544745 0;1 0.764450940710821 0;1 0.753505755876897 0;1 0.742560571042973 0;1 0.731615386209049 0;1 0.720670201375125 0;1 0.709725016541202 0;1 0.698779831707278 0;1 0.687834646873354 0;1 0.67688946203943 0;1 0.665944277205507 0;1 0.654999092371583 0;1 0.644053907537659 0;1 0.633108722703735 0;1 0.622163537869811 0;1 0.611218353035888 0;1 0.600273168201964 0;1 0.58932798336804 0;1 0.578382798534116 0;1 0.567437613700192 0;1 0.556492428866269 0;1 0.545547244032345 0;1 0.534602059198421 0;1 0.523656874364497 0;1 0.512711689530573 0;1 0.50176650469665 0;1 0.490821319862726 0;1 0.479876135028802 0;1 0.468930950194878 0;1 0.457985765360954 0;1 0.447040580527031 0;1 0.436095395693107 0;1 0.425150210859183 0;1 0.414205026025259 0;1 0.403259841191336 0;1 0.392314656357412 0;1 0.381369471523488 0;1 0.370424286689564 0;1 0.35947910185564 0;1 0.348533917021717 0;1 0.337588732187793 0;1 0.326643547353869 0;1 0.315698362519945 0;1 0.304753177686021 0;1 0.293807992852098 0;1 0.282862808018174 0;1 0.27191762318425 0;1 0.260972438350326 0;1 0.250027253516402 0;1 0.239082068682479 0;1 0.228136883848555 0;1 0.217191699014631 0;1 0.206246514180707 0;1 0.195301329346783 0;1 0.18435614451286 0;1 0.173410959678936 0;1 0.162465774845012 0;1 0.151520590011088 0;1 0.140575405177164 0;1 0.129630220343241 0;1 0.118685035509317 0;1 0.107739850675393 0;1 0.0967946658414692 0;1 0.0862782614443353 0;1 0.0706089693418639 0;1 0.0389084710450919 0;0.995823922717401 0.0113840500309187 0;0.975507474451548 0 0;0.943806976154776 0 0;0.912106477858004 0 0;0.880405979561232 0 0;0.84870548126446 0 0;0.817004982967688 0 0;0.785304484670916 0 0;0.753603986374144 0 0;0.721903488077372 0 0;0.6902029897806 0 0;0.658502491483828 0 0;0.626801993187056 0 0;0.595101494890285 0 0;0.563400996593513 0 0;0.531700498296741 0 0;0.500000000000003 0 0],...
tempcolornodeall=zeros(roi_num,size(Gradeintsall,2));
tempageall=zeros(roi_num,1);
tempcolornodeall(indnodes,:)=Gradeintsall;
tempageall(indnodes,1)=intersectionpointall;

indsig=find(allp<=pthref);
indnosig=find(allp>pthref);
ind1=find(allm==1);
ind2=find(allm==2);


indcall=[];
[xx numone]=sort(PerVonE(:,1),'descend');
for i=1:max(AALVonE(:,2))
indc=find(AALVonE(:,2)==i);
%indc1=intersect(indc,indnosig);
indc1=[];
indcraw2=intersect(intersect(indc,indsig),ind1);
indcraw3=intersect(intersect(indc,indsig),ind2);
[Vtemp Vdtemp]=sort(tempcolornodeall(indcraw2,1),'descend');
indc2= indcraw2(Vdtemp);
[Vtemp Vdtemp]=sort(tempcolornodeall(indcraw3,1),'descend');
indc3=indcraw3(Vdtemp);
indcall{i,1}=[indc1; indc3; indc2;];
indcall{i,2}=[length(indc1) length(indc3) length(indc2)];
end
indcalltogether=[];
for i=1:length(numone)
indcalltogether=[indcalltogether; indcall{numone(i),1}];
end
tempcolornodeplot=tempcolornode(indcalltogether,:);
tempcolornodeplot=tempcolornodeall(indcalltogether,:);
tempageallplot=tempageall(indcalltogether,:);

figure;imagesc(tempcolornodeplot.*tempcolornodeplot);caxis([0.000001 12]);hold on; 
%axis equal;
colormap(myslopecmap);
colorbar;
clp=0;clpall=[];
% for i=1:length(indcall)
%     clp=length(indcall{numone(i),1})+clp;
%     clpall(i)=clp;
%     plot([1:length(tempcolornodeplot)],repelem(clp,length(tempcolornodeplot)),'w','LineWidth',1.2);
%     drawnow
% end
plot(tempageallplot,[1:size(tempcolornodeplot,1)],'w','LineWidth',3.5);
hold off;
set(gca,'YTick',(clpall-5),'YTickLabel',[sysname(2,numone)]);
set(gca,'XTick',[Fitted_Subgroup],'XTickLabel',[0 0.5 1 1.5 2 2.5 3 3.5 4 5 6 8]);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',12,'LineWidth',0.1,'FontWeight','bold');
print(gcf,'-dtiff','-r300',[result_name 'SlopeVonEc_MA_FDRnn.tif']); close all;


indcall=[];
[xx numone]=sort(PerYeo(:,1),'descend');
for i=1:max(AALYeo(:,2))
indc=find(AALYeo(:,2)==i);
%indc1=intersect(indc,indnosig);
indc1=[];
indcraw2=intersect(intersect(indc,indsig),ind1);
indcraw3=intersect(intersect(indc,indsig),ind2);
[Vtemp Vdtemp]=sort(tempcolornodeall(indcraw2,1),'descend');
indc2= indcraw2(Vdtemp);
[Vtemp Vdtemp]=sort(tempcolornodeall(indcraw3,1),'descend');
indc3=indcraw3(Vdtemp);
indcall{i,1}=[indc1; indc3; indc2;];
indcall{i,2}=[length(indc1) length(indc3) length(indc2)];
end
indcalltogether=[];
for i=1:length(numone)
indcalltogether=[indcalltogether; indcall{numone(i),1}];
end
tempcolornodeplot=tempcolornode(indcalltogether,:);
tempcolornodeplot=tempcolornodeall(indcalltogether,:);
tempageallplot=tempageall(indcalltogether,:);

figure;imagesc(tempcolornodeplot);caxis([0.000001 5.5]);hold on; 
%axis equal;
colormap(myslopecmap);
colorbar;
clp=0;clpall=[];
% for i=1:length(indcall)
%     clp=length(indcall{numone(i),1})+clp;
%     clpall(i)=clp;
%     plot([1:length(tempcolornodeplot)],repelem(clp,length(tempcolornodeplot)),'w','LineWidth',1.2);
%     drawnow
% end
plot(tempageallplot,[1:size(tempcolornodeplot,1)],'w','LineWidth',1.2);
hold off;
set(gca,'YTick',(clpall-5),'YTickLabel',[sysname(3,numone)]);
set(gca,'XTick',[Fitted_Subgroup],'XTickLabel',[1 2  4 6 8]);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
print(gcf,'-dtiff','-r300',[result_name 'SlopeYeoc_MA_FDRnn.tif']); close all;



%VonE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VonE system slope stats 
pVonEplotall=[];
[xx numone]=sort(PerVonE(:,1),'descend');
% 
% for agei=1:length(Gradeintsall);
% %for agei=length(Gradeintsall);
% colornode(indnodes)=Gradeintsall(:,agei);
% colornodemean=zeros(500,1);colornodemediation=zeros(500,1);
% colornodemean(indnodes)=Mean_rapid_slope;
% colornodemediation(indnodes)=Median_rapid_slope;
% 
% 
% PerVonEslop=[];PerVonEmslop=[];PerVonEmdslop=[];
% for i=1:max(AALVonE(:,2))
%     tempm=allm(find(AALVonE(:,2)==i));  
%     tempp=allp(find(AALVonE(:,2)==i)); 
%     tempslope=colornode(find(AALVonE(:,2)==i));
%     tempmSig=tempslope(find(tempp<=pthref)); 
%     tempslopem=colornodemean(find(AALVonE(:,2)==i));
%     tempmSigm=tempslopem(find(tempp<=pthref)); 
%     tempslopemd=colornodemediation(find(AALVonE(:,2)==i));
%     tempmSigmd=tempslopemd(find(tempp<=pthref)); 
%     PerVonEslop{i}=tempmSig;
%     PerVonEmslop{i}=tempmSigm;
%     PerVonEmdslop{i}=tempmSigmd;
% end
% 
% 
% pVonE=[];QVon=[]; pVonEm=[];QVonm=[]; pVonEmd=[];QVonmd=[];
% for  it=1:7;
%     for jt=1:7;
%         if(it~=jt)
%     [pVonE(it,jt,agei), QVon(it,jt,agei)]= ranksum(PerVonEslop{1,it},PerVonEslop{1,jt});
%         [pVonEm(it,jt,agei), QVonm(it,jt,agei)]= ranksum(PerVonEmslop{1,it},PerVonEmslop{1,jt});
%             [pVonEmd(it,jt,agei), QVonmd(it,jt,agei)]= ranksum(PerVonEmdslop{1,it},PerVonEmdslop{1,jt});
%         else
%           pVonE(it,jt,agei)=999;  QVon(it,jt,agei)=999;
%                     pVonEm(it,jt,agei)=999;  QVonm(it,jt,agei)=999;
%                               pVonEmd(it,jt,agei)=999;  QVonmd(it,jt,agei)=999;
%         end
%     end
% end
% 
% 
% pVonEnew=pVonE(numone,numone,agei);
% pVonEplot=-log(pVonEnew);
% pVonEplot(find(pVonEnew==999))=0;
% pVonEmnew=pVonEm(numone,numone,agei);
% pVonEmplot=-log(pVonEmnew);
% pVonEmplot(find(pVonEmnew==999))=0;
% pVonEmdnew=pVonEmd(numone,numone,agei);
% pVonEmdplot=-log(pVonEmdnew);
% pVonEmdplot(find(pVonEmdnew==999))=0;
% pall=[];
% pall=pVonEnew(intersect((find(triu(pVonEnew)>0)),(find(triu(pVonEnew)<1))));
% [pmmfdr1 pmmfdr2]=gretna_FDR(pall,0.05);
% if(isempty(pmmfdr1))
%     pVonEplot=zeros(size(numone,1),size(numone,1));
% else
% pVonEplot(find(pVonEnew>=(pmmfdr1)))=0;
% end
% pVonEplotall(:,:,agei)=pVonEplot;
% pVonEmplotall(:,:,agei)=pVonEmplot;
% pVonEmdplotall(:,:,agei)=pVonEmdplot;
% pVonEmnewall(:,:,agei)=pVonEmnew;
% pVonEmdnewall(:,:,agei)=pVonEmdnew;
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %VonE system slope stats  Pmap plot
% % 创建一个示例的 3D 密度分布矩阵，这里用随机数据代替
% density_matrix = pVonEplotall(:,:,1:end);
% density_matrixnew=[];
% for i=1:length(density_matrix)
%     density_matrixnew(:,:,i)=density_matrix(:,:,(length(density_matrix)-i)+1);
% end
% % 创建一个新的 figure
% figure;
% 
% % 创建切片图像
% h = slice(density_matrixnew, [], [], 1:size(density_matrixnew, 3));
% 
% % 设置绘图属性和标签
% zlabel('Age');
% set(gca,'XTick',[1:1:length(sysname(2,numone))],'XTickLabel',[sysname(2,numone)]);
% set(gca,'YTick',[1:1:length(sysname(2,numone))],'YTickLabel',([sysname(2,numone)]));
% set(gca,'ZTick',flip([length(density_matrixnew)-Fitted_Subgroup]),'ZTickLabel',[8 6 5 4 3.5 3 2.5 2 1.5 1 0.5 0]);
% title('P matrix of slopes');
% % 设置视角
% view(3); % 3D 视角
% 
% % 移除边缘线
% for i = 1:numel(h)
%     set(h(i), 'EdgeColor', 'none');
% end
% 
% set(gcf,'Color',[1 1 1]) ;
% set(gcf,'position', [1 70 970 800]);
% set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',12,'LineWidth',0.1,'FontWeight','bold');
% print(gcf,'-dtiff','-r300',[result_namethre 'Slopelinematrix_VonEc.tif']); close all;


%pmap 圆形 平均速率

% for agei=1:length(Gradeintsall);
for agei=1;
colornode(indnodes)=Gradeintsall(:,agei);
colornodemean=zeros(500,1);colornodemediation=zeros(500,1);
colornodemean(indnodes)=Mean_rapid_slope;
colornodemediation(indnodes)=Median_rapid_slope;


PerVonEslop=[];PerVonEmslop=[];PerVonEmdslop=[];
for i=1:max(AALVonE(:,2))
    tempm=allm(find(AALVonE(:,2)==i));  
    tempp=allp(find(AALVonE(:,2)==i)); 
    tempslope=colornode(find(AALVonE(:,2)==i));
    tempmSig=tempslope(find(tempp<=pthref)); 
    tempslopem=colornodemean(find(AALVonE(:,2)==i));
    tempmSigm=tempslopem(find(tempp<=pthref)); 
    tempslopemd=colornodemediation(find(AALVonE(:,2)==i));
    tempmSigmd=tempslopemd(find(tempp<=pthref)); 
    PerVonEslop{i}=tempmSig;
    PerVonEmslop{i}=tempmSigm;
    PerVonEmdslop{i}=tempmSigmd;
end


pVonE=[];QVon=[]; pVonEm=[];QVonm=[]; pVonEmd=[];QVonmd=[];
for  it=1:7;
    for jt=1:7;
        if(it~=jt)
    [pVonE(it,jt,agei), QVon(it,jt,agei)]= ranksum(PerVonEslop{1,it},PerVonEslop{1,jt});
        [pVonEm(it,jt,agei), QVonm(it,jt,agei)]= ranksum(PerVonEmslop{1,it},PerVonEmslop{1,jt});
            [pVonEmd(it,jt,agei), QVonmd(it,jt,agei)]= ranksum(PerVonEmdslop{1,it},PerVonEmdslop{1,jt});
        else
          pVonE(it,jt,agei)=999;  QVon(it,jt,agei)=999;
                    pVonEm(it,jt,agei)=999;  QVonm(it,jt,agei)=999;
                              pVonEmd(it,jt,agei)=999;  QVonmd(it,jt,agei)=999;
        end
    end
end


pVonEnew=pVonE(numone,numone,agei);
pVonEplot=-log(pVonEnew);
pVonEplot(find(pVonEnew==999))=0;
pVonEmnew=pVonEm(numone,numone,agei);
pVonEmplot=-log(pVonEmnew);
pVonEmplot(find(pVonEmnew==999))=0;
pVonEmdnew=pVonEmd(numone,numone,agei);
pVonEmdplot=-log(pVonEmdnew);
pVonEmdplot(find(pVonEmdnew==999))=0;
pall=[];
pall=pVonEmnew(intersect((find(triu(pVonEmnew)>0)),(find(triu(pVonEmnew)<1))));
[pmmfdr1 pmmfdr2]=gretna_FDR(pall,0.05);
pthrenode=0.05;
if(isempty(pmmfdr1))
    pVonEplot=zeros(size(numone,1),size(numone,1));
else
pVonEplot(find(pVonEnew>=(pthrenode)))=0;
pVonEmplot(find(pVonEmnew>=(pthrenode)))=0;
pVonEmdplot(find(pVonEmdnew>=(pthrenode)))=0;
end
pVonEplotall(:,:,agei)=pVonEplot;
pVonEmplotall(:,:,agei)=pVonEmplot;
pVonEmdplotall(:,:,agei)=pVonEmdplot;
pVonEmnewall(:,:,agei)=pVonEmnew;
pVonEmdnewall(:,:,agei)=pVonEmdnew;
end

% 
% pVonEmnew=pVonEmnewall(:,:,end);
% pallm=pVonEmnew(intersect((find(triu(pVonEmnew)>0)),(find(triu(pVonEmnew)<1))));
% [pmmmfdr1 pmmmfdr2]=gretna_FDR(pallm,0.05);
% pmmmfdr1=0.01;
% pVonEmplot=pVonEmplotall(:,:,end);
% pVonEmplot(find(pVonEmnew>pmmmfdr1))=0;
% 
% pVonEmdnew=pVonEmdnewall(:,:,end);
% pallmd=pVonEmdnew(intersect((find(triu(pVonEmdnew)>0)),(find(triu(pVonEmdnew)<1))));
% [pmmmdfdr1 pmmmdfdr2]=gretna_FDR(pallmd,0.05);
% pmmmdfdr1=0.01;
% pVonEmdplot=pVonEmdplotall(:,:,end);
% pVonEmdplot(find(pVonEmdnew>pmmmdfdr1))=0;

% 创建一个绘图函数，只在圆内填充伪彩色的颜色
figure;
data=(pVonEmplot);
hold on;
[row, col] = size(data);
[x, y] = meshgrid(1:col, 1:row);
radius = 0.4; % 圆的半径

% 创建一个颜色映射
cmap = nclCM(365);
cmap(1,:)=[0.89 0.89 0.89];
for i = 1:row
    for j = 1:col
        % 计算圆心的坐标
        cx = j;
        cy = i;
        
        % 计算颜色索引，限制在1到256之间
        color_idx = round(data(i, j)./max(max(data)) * (length(cmap)-1))+1;
        
        % 从颜色映射中获取颜色
        color = cmap(color_idx, :);
        
        rectangle('Position', [j - radius, i - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', color, 'EdgeColor', 'none');
    end
end
axis equal;
colormap(cmap);
colorbar;
axis ij %// use if you want the origin at the top, like with imagesc
set(gca,'XTick',[1:1:length(pVonEmplot)],'XTickLabel',[sysname(2,numone)]);
set(gca,'YTick',[1:1:length(pVonEmplot)],'YTickLabel',[sysname(2,numone)]);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
%axis off
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
hold off;
print(gcf,'-dtiff','-r300',[result_namethre 'PerVonEP_MeanSlope_Circlemap.tif']); close all;

figure;
data=(pVonEmdplot);
hold on;
[row, col] = size(data);
[x, y] = meshgrid(1:col, 1:row);
radius = 0.4; % 圆的半径

% 创建一个颜色映射
cmap = nclCM(365);
cmap(1,:)=[0.89 0.89 0.89];
for i = 1:row
    for j = 1:col
        % 计算圆心的坐标
        cx = j;
        cy = i;
        
        % 计算颜色索引，限制在1到256之间
        color_idx = round(data(i, j)./max(max(data)) * (length(cmap)-1))+1;
        
        % 从颜色映射中获取颜色
        color = cmap(color_idx, :);
        
        rectangle('Position', [j - radius, i - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', color, 'EdgeColor', 'none');
    end
end
axis equal;
colormap(cmap);
colorbar;
axis ij %// use if you want the origin at the top, like with imagesc
set(gca,'XTick',[1:1:length(pVonEmdplot)],'XTickLabel',[sysname(2,numone)]);
set(gca,'YTick',[1:1:length(pVonEmdplot)],'YTickLabel',[sysname(2,numone)]);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
%axis off
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',20,'LineWidth',0.1,'FontWeight','bold');
hold off;
print(gcf,'-dtiff','-r300',[result_namethre 'PerVonEP_MedianSlope_Circlemap.tif']); close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VonE system slope bar stats plot

%%von
meanall=[];
for ii=1:size(PerVonEmslop,2)
meanalltemp=mean(PerVonEmslop{1,ii})
meanall=[meanall meanalltemp];
end
%[xx numone]=sort((meanall));
%numone=VonEnumonefinal;
% PerVonEslopnew=[];
% for ii=1:size(PerVonEmslop,2)
% PerVonEslopnew{1,ii}=PerVonEmslop{1,numone(ii)};
% end
x=[];g=[];
for ii=1:size(PerVonEmslop,2)
x = [x; PerVonEmslop{1,ii}];
gt = repmat({[sysname{2,ii}]},length(PerVonEmslop{1,ii}),1);
g=[g; gt];
end

figure
%boxplot(x,g,'Notch','on')
violinplot(x,g,'GroupOrder',sysname(2,numone),'ShowData',false,'ShowMean',true,'ShowNotches',true);
%barplot(PerHer4new(:,1)+PerHer4new(:,2)+PerHer4new(:,3),'FaceColor',[70/255 171/255,215/255]);hold on;barh(PerHer4new(:,2)+PerHer4new(:,3),'FaceColor',[237/255 176/255 33/255]);hold on;barh(PerHer4new(:,3),'FaceColor',[152/255 152/255 152/255]);
set(gca,'XTick',[1:1:length(PerVonEmslop)],'XTickLabel',[sysname(2,numone)]);
ylabel('Velocity');
set(gcf,'position', [4         20       1000        900]);
set(gca,'position',[0.1,0.22,0.8,0.75]);
set(gcf, 'PaperPositionMode', 'auto');set(gca, 'fontsize',25,'LineWidth',0.1,'FontWeight','bold');
print(gcf,'-dtiff','-r300',[result_namethre 'VonE_meanSlope_voilin.tif']); close all;


%%%%%%%%%%%%%%%%%%%%%
% %SVR
% % load([Pathnew(1:end-20),'\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);
% % nodal_eff=Ne';

% nodal_eff=Nodal_metric';
% finalYCorr=mean(nodal_eff(:,find((Age>7.5)&(Age<=9))),2);
load([Pathnewplot 'Forprediction_finalYCorr.mat'])
dirsvr=(Pathnewplot);
regionalcorrsvr(indnodes,Fitted_Subgroup,Age2,roi_num,Gradeintsall,finalYCorr,typeofstat,result_namethre,'SVR_Ageat8ys',[],dirsvr);




end


