%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% brain hubs
% set working paths
statdir=[Pathnew,'Forgithub\Code'];
savedir=[Pathnew,'Forgithub'];
Pathnewplot=[Pathnew,'Forgithub'];
cd(Pathnewplot)
load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);
plotornot=1;


%% sliding windows
[sort_age sort_age_ind]=sort(Age_gender(:,1));
sub_ind=[];
validgn=10;
parc=1-(1./validgn);
round(114*parc);


Length_slide_windows=15;
slideG=[];meanageslide=[];slideG_order=[];
for i=1:length(Age_gender)+1-Length_slide_windows;
slideG_order(:,i)=[1:Length_slide_windows]+i-1;
meanageslide(i)=mean(sort_age(slideG_order(:,i)))./365;
slideG(:,i)=sort_age_ind(slideG_order(:,i));
end


Sv1=13;Sv2=12;
Nname='NgE';
Nodalmetrics=Ne;
result_name='NgE';

hdrn = spm_vol([Pathnew,'Forgithub\Data\500parcellation.nii']);
Atlas=spm_read_vols(hdrn);
hdrn.dt=[16 1];

hdre=hdrn;hdrep=hdre;
ShowG=[0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];ShowGind=[];
%ShowG=[0.5 1.0 1.5 2 2.5 3 3.5];ShowGind=[];
for i=1:length(ShowG)
[temp ShowGind(i)]=min(abs(meanageslide-ShowG(i)));
end
ShowGindnew = unique(ShowGind);
meanageslide(ShowGindnew);

XTicklabelstr=[];
for i=1:length(ShowGindnew)
XTicklabelstr=[XTicklabelstr;num2str(meanageslide(ShowGindnew(i)),'%4.2f')];
end
slideGG=ShowGindnew;
minmaxvalue=[];minmaxvalueforplot=[];
for ii=1:length(ShowGindnew);
Atlaszeros2=zeros(size(Atlas));Atlaszeros3=zeros(size(Atlas));Atlaszeros4=zeros(size(Atlas));Atlaszeros2forplot=zeros(size(Atlas));
for jj=1:length(unique(Atlas(find(Atlas~=0))));
Atlaszeros2(find(Atlas==jj))=mean(Nodalmetrics(slideG(:,ShowGindnew(ii)),jj));
end
minmaxvalue(ii,:)=[min(Atlaszeros2(find(Atlaszeros2~=0))) max(Atlaszeros2(find(Atlaszeros2~=0)))];
hdre.fname=[result_name, '_SlideGG_', XTicklabelstr(ii,:),'.nii'];
spm_write_vol(hdre,Atlaszeros2);

Atlaszeros2forplot(find(Atlaszeros2~=0))=(Atlaszeros2(find(Atlaszeros2~=0)));
minmaxvalueforplot(ii,:)=[min(Atlaszeros2forplot(find(Atlaszeros2forplot~=0))) max(Atlaszeros2forplot(find(Atlaszeros2forplot~=0)))];
hdrep.fname=[result_name '_SlideGGforplot_', XTicklabelstr(ii,:),'.nii'];
spm_write_vol(hdrep,Atlaszeros2forplot);
uniqevalue=unique(Atlaszeros2(find(Atlaszeros2~=0)));[sortv sorto]=sort(uniqevalue);

for jj=1:length(unique(Atlas(find(Atlas~=0))));
Atlaszeros3(find(Atlaszeros2==sortv(jj)))=sorto(jj);
end

hdrn.fname=[result_name 'Rank_SlideGG' num2str(meanageslide(slideGG(ii))) '.nii'];
spm_write_vol(hdrn,Atlaszeros3);
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',hdrn.fname,[Pathnew,'Forgithub\Code\opt\hubrankcolor_opt.mat'],[result_name 'Eff_SlideGG' num2str(meanageslide(slideGG(ii))) '.tif']);
close all
end


%% clustering
Nodalmetrics=Ne;
slideGG=[1:1:length(meanageslide)];
minmaxvalueall=[];meaneffslide=[];meaneffrankslide=[];
for ii=1:length(slideGG);
for jj=1:length(unique(Atlas(find(Atlas~=0))));
meaneffslide(jj,ii)=mean(Nodalmetrics(slideG(:,slideGG(ii)),jj));
end
[sortvr sortor]=sort(meaneffslide(:,ii));
meaneffrankslide(:,ii)=sortor;
minmaxvalueall(ii,:)=[min(meaneffslide(:,ii)) max(meaneffslide(:,ii))];
end


per=30;
nodalclusternewm=meaneffslide;
[pc,score,latent,tsquare] = pca(nodalclusternewm');
perpca=cumsum(latent)./sum(latent);
tran=pc(:,1:per);
nodalEg_All_Nodenorm_afterpca= bsxfun(@minus,nodalclusternewm',mean(nodalclusternewm',1));
nodalEg_All_Nodenorm_afterpca= nodalEg_All_Nodenorm_afterpca*tran;
[maxpar nc]=nodalclusternew(nodalEg_All_Nodenorm_afterpca','cosine','sample',5,100,'BarclusterEffAge.tif');
Agepeakclusterindex=find(diff(maxpar)~=0);
ShowG=[0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];ShowGind=[];

for i=1:length(ShowG)
[temp, ShowGind(i)]=min(abs(meanageslide-ShowG(i)))
end
ShowGindnew = unique(ShowGind);
meanageslide(ShowGindnew);
figure;plot(maxpar,'v','LineWidth',1.5,'MarkerFaceColor',[220/255 11/255 50/255],'MarkerSize',8,'MarkerEdgeColor',[180/255 11/255 50/255]);
hold on;
XTicklabelstr=[];
for i=1:length(ShowGindnew)
XTicklabelstr=[XTicklabelstr;num2str(meanageslide(ShowGindnew(i)),'%4.1f')];
end
set(gca,'XTick',ShowGindnew','XTicklabel',XTicklabelstr);
set(gca,'YTick',[0 1 2],'YTicklabel',['0'; '1'; '2']);
set(gcf,'position', [12         12      1300        600]);
set(gca,'position',[0.1,0.1,0.85,0.85]);
set(gcf, 'PaperPositionMode', 'auto')  
set(gca, 'fontsize', 18,'LineWidth',2.5,'FontWeight','bold');
ylabel('Cluster');

print(gcf,'-dtiff','-r300',['ClusterEffAge',num2str(per),'.tif']);



figure;
Matrixcorr=corr(nodalclusternewm-repmat(mean(nodalclusternewm,2),1,100));
imagesc(Matrixcorr(1:2:end,1:2:end));colormap(hot);hold on;
set(gca,'XTick',floor(ShowGind(1:end-1)./2),'XTicklabel',['0.5 y'; '  1 y'; '1.5 y'; '  2 y'; '2.5 y'; '  3 y';'  4 y'; '  5 y'; '  6 y';'  7 y';]);
xtickangle(90)
set(gca,'YTick',floor(ShowGind(1:end-1)./2),'YTicklabel',['0.5 y'; '  1 y'; '1.5 y'; '  2 y'; '2.5 y'; '  3 y';'  4 y'; '  5 y'; '  6 y';'  7 y';]);
set(gcf,'position', [12     12    900       850]);
set(gca,'position',[0.15 0.15,0.8,0.8]);
set(gcf, 'PaperPositionMode', 'auto')  
set(gca, 'fontsize', 30,'LineWidth',2.5,'FontWeight','bold');
print(gcf,'-dtiff','-r300','ClusterEffAgeMatrixCos.tif');

minmaxvaluecluster=[];
for ii=1:max(maxpar);
Atlaszeros2=zeros(size(Atlas));Atlaszeros3=zeros(size(Atlas));Atlaszeros4=zeros(size(Atlas));
for jj=1:length(unique(Atlas(find(Atlas~=0))));
Atlaszeros2(find(Atlas==jj))=mean(Nodalmetrics(sort_age_ind(find(maxpar==ii)),jj));
end
minmaxvaluecluster(ii,:)=[min(Atlaszeros2(find(Atlaszeros2~=0))) max(Atlaszeros2(find(Atlaszeros2~=0)))];
hdre.fname=[result_name '_SlideCluster' num2str(ii) '.nii'];
spm_write_vol(hdre,Atlaszeros2);
uniqevalue=unique(Atlaszeros2(find(Atlaszeros2~=0)));[sortv sorto]=sort(uniqevalue);
for jj=1:length(unique(Atlas(find(Atlas~=0))));
Atlaszeros3(find(Atlaszeros2==sortv(jj)))=sorto(jj);
end
hdrn.fname=[result_name 'Rank__SlideCluster' num2str(ii) '.nii'];
spm_write_vol(hdrn,Atlaszeros3);
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',hdrn.fname,[Pathnew,'Forgithub\Code\opt\hubrankcolor_opt.mat'],[result_name 'Rank__SlideCluster' num2str(ii)  '.tif']);
close all
end

hdrnr1=spm_vol([result_name 'Rank__SlideCluster1.nii']);
atlasnr1=spm_read_vols(hdrnr1);
hdrnr2=spm_vol([result_name 'Rank__SlideCluster2.nii']);
atlasnr2=spm_read_vols(hdrnr2);
atlasnr1(find(atlasnr1>424))=1000;
atlasnr1(find(atlasnr1<=424))=1;
atlasnr2(find(atlasnr2>424))=2000;
atlasnr2(find(atlasnr2<=424))=1;
atlasnr3=atlasnr2-atlasnr1;
hdrnr2.fname=[result_name 'Rank__SlideClusterall.nii'];
spm_write_vol(hdrnr2,atlasnr3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% brain modules
clear; clc
load([Pathnew,'Forgithub\Data\Network\Matrixallpar500FNsortage.mat'); 
intralayer_resolution = 1; % intralayer resolution for multilayer modularity
interlayer_resolution = 1; % interlayer resolution for multilayer modularity
N_nodes = 500;
for i = 1:100
    AA = Matrixallpar500FNsortage;
    [B,mm] = multiord(AA,intralayer_resolution,interlayer_resolution); % B=supra-adjacency matrix
    PP = @(S) postprocess_ordinal_multilayer(S,114); % multilayer modularity set-up
    [S,Q1,mod_iter_tmp] = iterated_genlouvain(B,10000,0,1,'moverandw',[],PP); % S = modularity assignments; Q1 = modularity value; mod_iter_tmp = iterations
    S = reshape(S,N_nodes,114); % 2D representation of modularity assignments
    modularity_assignments{i} = S;
    number_of_communities(i,:) = sum(diff(sort(S)) ~= 0) +1; % total number of community in network
end

load([Pathnew,'Forgithub\Data\Network\Multilayer_modules_norm_100.mat'])

%% slidewindows
Agesubind=[];
[sort_age sort_age_ind]=sort(Age_gender(:,1));
sub_ind=[];
validgn=10;
parc=1-(1./validgn);
round(114*parc);
% 
% 
Length_slide_windows=15;
slideG=[];meanageslide=[];slideG_order=[];
for i=1:length(brainsize)+1-Length_slide_windows;
slideG_order(:,i)=[1:Length_slide_windows]+i-1;
meanageslide(i)=mean(sort_age(slideG_order(:,i)))./365;
slideG(:,i)=sort_age_ind(slideG_order(:,i));
end

Minnz=[];Mimzr=[];

for i=1:8;for j=1:8;Minnz(i,j)=Allznmi{i,j}(1);end;end;
for i=1:8;for j=1:8;Mimzr(i,j)=Allznmi{i,j}(4);end;end;
figure
imagesc(Minnz)
figure
imagesc(Mimzr)


for i=1:8;for j=1:8;Meanz(i,j)=Allznmi{i,j}(3);end;end;
for i=1:8;for j=1:8;Meanzr(i,j)=Allznmi{i,j}(6);end;end;
Meanzc=Meanzr-Meanz;
figure
imagesc(Meanz)
figure
imagesc(Meanzr)
figure
imagesc(Meanzc)


meanncom=[];
for i=1:8
    for j=1:8
ncom=squeeze(number_of_communitiesall(i,j,:,:));

meanncom{i,j}=mean(ncom,1);
    end
end


tt=2;tj=1;
AllQf=Q1allall(tt,tj,:);
[ma md]=max(AllQf)
ncf=squeeze(number_of_communitiesall(tt,tj,md,:))
plot(ncf)

ncfall=squeeze(number_of_communitiesall(tt,tj,:,:));
meanncfall=mean(ncfall);
meanncfallnoover=meanncfall(1:15:100);
meanageslidenew=meanageslide(1:15:100);
[r p]=corr(meanageslidenew',meanncfallnoover')
plot(meanncfallnoover)

[r p]=corr(meanageslide(1:1:100)',ncf(1:1:100))


modularity_assignmentsselect=[];
for i=1:100;
modularity_assignmentsselect{i}=modularity_assignmentsall{tt,tj,i};
end
k=md
ShowG=[0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];ShowGind=[];
for i=1:length(ShowG)
[temp, ShowGind(i)]=min(abs(meanageslide-ShowG(i)))
end
ShowGindnew = unique(ShowGind);
meanageslide(ShowGindnew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot modules
numofroi=500;
hdrn = spm_vol([Pathnew, 'Forgithub\Data\500parcellation.nii']);
Atlas=spm_read_vols(hdrn);
for iss =unique(ShowGind)
result_name=['2023NMModallSlidewdGW15_',num2str(tt),'_',num2str(tj),'_',num2str(k),'_',num2str(iss,'%03d'),''];
AtlasPeak=zeros(size(Atlas));hdrn.dt=[16 1];
if(iss==1)
modularity_assig1=modularity_assignmentsselect{1,k}(:,iss);
%modularity_assig1=M.Ci;
for i=1:numofroi;
    AtlasPeak(find(Atlas==i))=modularity_assig1(i);
end

hdrn.fname=[result_name '.nii'];
spm_write_vol(hdrn,AtlasPeak);
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',hdrn.fname,[Pathnew,'Forgithub\Code\opt\module_opt.mat'],[ result_name '.tif']);

else
modularity_assigleft=modularity_assignmentsselect{1,k}(:,iss);
for i=1:numofroi;
    AtlasPeak(find(Atlas==i))=modularity_assigtemp(i);
end

for i=1:numofroi;
    AtlasPeak(find(Atlas==i))=modularity_assigleft(i);
end

hdrn.fname=[result_name '.nii'];
spm_write_vol(hdrn,AtlasPeak);
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',hdrn.fname,[Pathnew,'Forgithub\Code\opt\module_opt.mat'],[result_name '.tif']);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot module flexibility
IndLabel=modularity_assignmentsselect{1,k};
[V_wei, V_vec, Deg_J] = scaled_inclusivity_wei(IndLabel',3);


hdrn = spm_vol([Pathnew,'Forgithub\Data\500parcellation.nii']);
Atlas=spm_read_vols(hdrn);
AtlasPeak=zeros(size(Atlas));
hdrn.dt=[16 1];
for i=1:numofroi;
    AtlasPeak(find(Atlas==i))=V_wei(i);
end
result_name='MV';
hdrn.fname=[result_name '.nii'];
spm_write_vol(hdrn,AtlasPeak);
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',hdrn.fname,[Pathnew,'Forgithub\Code\opt\mv_opt.mat'],[result_name '.tif']);


modularity_assignmentnew=modularity_assignmentsselect;
selectmod=modularity_assignmentsselect{k};

IndLabelFinalkallmatch=zeros(500,size(modularity_assignmentnew,2));
firstnodemod=selectmod(:,1);
for v=2:size(modularity_assignmentnew,2)
tempnodemod=selectmod(:,v);
[IndLabelFinalkallmatch(:,v)] = BestMapping(firstnodemod,tempnodemod);
end
IndLabelFinalkallmatch(:,1)=firstnodemod;

[V_weiO, V_vecO, Deg_JO] = scaled_inclusivity_wei(IndLabelFinalkallmatch',3);





