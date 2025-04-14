% set working paths
statdir=[Pathnew,'Forgithub\Code'];
savedir=[Pathnew,'Forgithub'];
Pathnewplot=[Pathnew,'Forgithub'];
cd(Pathnewplot)
load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);
plotornot=1;
load([Pathnew,'Forgithub\Data\Network\Matrixallpar500FN.mat']);

load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);

md4plot=indminNodalStrmdall{1}{1,1}{1,1};
    tbl=md4plot.Variables;%
    Age=tbl.Age;
nodal_eff=Ne';

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
    Ypoint=zeros(1,length(Index2))+((8/9).*(max(Index2)-min(Index2)))+min(Index2);
    Fp=(Ypoint)-Index2;
    Fp(find(Fp>0))=1;Fp(find(Fp<0))=-1;
    dffFp=diff(Fp);
    intersectionpointall(i)=find(dffFp~=0);
    intersectionAge(i)=Age2(intersectionpointall(i));
    Mean_rapid_slope(i)=mean(Gradeintsall(i,1:intersectionpointall(i)));
    Median_rapid_slope(i)=mean(Gradeintsall(i,1:intersectionpointall(i)));
    Index2;
end



NumPerm=10000; 
%% real slope
nodal_metricssig=zscore(Mean_rapid_slope(indnodes));
indnodesplot=indnodes;
corrstatsreal=[];
%The following images of brain axis should be warped into a common standard
%template by neuronmaps. 
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\Cortical.Thickness.dscalar.nii';
corrstatsreal{1}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','Thickness','GAM_Zs_Slope_Thickness')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\Evolution.Expansion.dscalar.nii';
corrstatsreal{2}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','Evolution','GAM_Zs_Slope_Expansion')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\CBF.dscalar.nii';
corrstatsreal{3}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','CBF','GAM_Zs_Slope_CBF')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\BigBrain.Histology.dscalar.nii';
corrstatsreal{4}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','Bigbrain Histology','GAM_Zs_Slope_BigbrainHistology')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\AllometricScaling.PNC20mm.dscalar.nii';
corrstatsreal{5}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','AllometricScaling','GAM_Zs_Slope_AllometricScaling')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\PC1.AHBA.dscalar.nii';
corrstatsreal{6}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','PC1.AHBA','GAM_Zs_Slope_PC1.AHBA')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\PC1.Neurosynth.dscalar.nii';
corrstatsreal{7}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','PC1.Neurosynth','GAM_Zs_Slope_PC1.Neurosynth')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\PET.AG.dscalar.nii';
corrstatsreal{8}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','AG','GAM_Zs_Slope_AG')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\G1.fMRI.dscalar.nii';
corrstatsreal{9}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','G1.fMRI','GAM_Zs_Slope_G1.fMRI')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\T1T2ratio.dscalar.nii';
corrstatsreal{10}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','T1T2ratio','GAM_Zs_Slope_T1T2ratio')
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\SensorimotorAssociation_Axis.dscalar.nii';
corrstatsreal{11}=surf_sig_roi_corr(nodal_metricssig,indnodesplot,scalarfile,'Nodal slope','SensorimotorAssociation_Axis','GAM_Zs_Slope_SensorimotorAssociation_Axis')

all_real_r=[];
for jj=1:11
   all_real_r(jj)= corrstatsreal{jj}{1,1};
end



nodal_metrics = Mean_rapid_slope;
lhgii=[Pathnew,'Forgithub\Data\Atlas\Par_500_fslr-lh.gii'];
rhgii=[Pathnew,'Forgithub\Data\Atlas\Par_500_fslr-rh.gii'];

LV=gifti(lhgii); % FS_LR 32k
RV=gifti(rhgii);
NumHemi=size(LV.cdata, 1);
AAL500LR=[LV.cdata; RV.cdata];
Perm_Results = zeros(500,NumPerm);

%% random test
corrstats1=[];
for ii=1:NumPerm
    if(mod(ii, 100)==0)
    disp(ii)
    end
    OutFile=['...\results\Rand_' sprintf('%.5d', ii)];
    %Random pars should be priorly generated by allnode_SpinTest.m or surrogate_spatial_resample.py
    load(OutFile);
    Rand_Label = [LNewLabel;RNewLabel];
    for node=1:500;
        index = find(AAL500LR==node);
        perm_roi = mode(Rand_Label(index));
        if perm_roi == 0
            Perm_Results(node,ii)=0;
            continue;
        end
        Perm_Results(node,ii) =nodal_metrics(perm_roi);
    end
nodal_metricssig_per=Perm_Results(indnodes,ii)';
indnodesplot=indnodes;

%The following images of brain axis should be warped into a common standard
%template by neuronmaps. 
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\Cortical.Thickness.dscalar.nii';
corrstats1{1,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','Thickness','Ted_Thickness',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\Evolution.Expansion.dscalar.nii';
corrstats1{2,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','Evolution','Ted_Expansion',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\CBF.dscalar.nii';
corrstats1{3,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','CBF','Ted_CBF',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\BigBrain.Histology.dscalar.nii';
corrstats1{4,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','Bigbrain Histology','Ted_BigbrainHistology',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\AllometricScaling.PNC20mm.dscalar.nii';
corrstats1{5,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','AllometricScaling','AllometricScaling',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\PC1.AHBA.dscalar.nii';
corrstats1{6,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','PC1.AHBA','PC1.AHBA',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\PC1.Neurosynth.dscalar.nii';
corrstats1{7,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','PC1.Neurosynth','PC1.Neurosynth',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\PET.AG.dscalar.nii';
corrstats1{8,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','AG','AG',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\G1.fMRI.dscalar.nii';
corrstats1{9,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','G1.fMRI','G1.fMRI',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\T1T2ratio.dscalar.nii';
corrstats1{10,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','T1T2ratio','T1T2ratio',Pathnew);
scalarfile='...\S-A_ArchetypalAxis-main\FSLRVertex\SensorimotorAssociation_Axis.dscalar.nii';
corrstats1{11,ii}=surf_sig_roi_corr(nodal_metricssig_per,indnodesplot,scalarfile,'Nodal slope','SensorimotorAssociation_Axis','SensorimotorAssociation_Axis',Pathnew);
end

allper_r1=[];
for ii=1:NumPerm
for jj=1:11
   allper_r1(jj,ii)= corrstats1{jj,ii}{1,1};
end
end

numnosig1=[];
for jj=1:11
    if(all_real_r(jj)>0)
   numnosig1(jj)=length(find(allper_r1(jj,:)>all_real_r(jj)));
    elseif(all_real_r(jj)<0)
   numnosig1(jj)=length(find(allper_r1(jj,:)<all_real_r(jj)));
    end
end

