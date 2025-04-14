function corrstats=surf_sig_roi_corr(Mean_rapid_slope,indnodes,cifti_scalar,Xlabels,Ylablels,namef,Pathnew)
lhgii=[Pathnew,'\Forgithub\Data\Atlas\Par_500_fslr-lh.gii'];
rhgii=[Pathnew,'\Forgithub\Data\Atlas\Par_500_fslr-rh.gii'];
%cifti_scalar=['D:\Documents\OneDrive - bnu.edu.cn\My_publications\0-8years\atlases_Oct2023\atlases_Jan1220\ted\S-A_ArchetypalAxis-main\FSLRVertex\Cortical.Thickness.dscalar.nii'];
thissurfl = gifti(lhgii);
labelssurfl= thissurfl.cdata;
labelssurflqni = unique(labelssurfl(find(labelssurfl>0)));
thissurfr = gifti(rhgii);
labelssurfr= thissurfr.cdata;
labelssurfrqni = unique(labelssurfr(find(labelssurfr>0)));

dsc_hemi=cifti_read(cifti_scalar);
leftscalartmp = cifti_struct_dense_extract_surface_data(dsc_hemi, 'CORTEX_LEFT');
rightscalartmp = cifti_struct_dense_extract_surface_data(dsc_hemi, 'CORTEX_RIGHT');

leftscalar=leftscalartmp(:,1);
rightscalar=rightscalartmp(:,1);

zerosindl=union(find(labelssurfl==0),find(leftscalar==0));
labelssurfl(zerosindl)=[];
leftscalar(zerosindl)=[];

zerosindr=union(find(labelssurfr==0),find(rightscalar==0));
labelssurfr(zerosindr)=[];
rightscalar(zerosindr)=[];
corrstats = corrstat(labelssurfl,leftscalar,[]);
%corrplot(labelssurfl,corrstats,'Mean roi','ddd',[],[],[],[],[],['new.tif'])

LROI_mean_scalar=[];
uqn_roil=unique(labelssurfl);
n_roi=length(uqn_roil);
for i=1:n_roi
    leftsca_n=leftscalar(find(labelssurfl==uqn_roil(i)));
    LROI_mean_scalar(i)=nanmean(leftsca_n(:));
    %LROI_mean_scalar(i)=mean(leftsca_n);
    
end

RROI_mean_scalar=[];
uqn_roir=unique(labelssurfr);
n_roi=length(uqn_roir);
for i=1:n_roi
    rightsca_n=rightscalar(find(labelssurfr==uqn_roir(i)));
    RROI_mean_scalar(i)=nanmean(rightsca_n(:));
   % RROI_mean_scalar(i)=mean(rightsca_n);
end

AllROI_m_sca = [LROI_mean_scalar RROI_mean_scalar];
All_label=[uqn_roil; uqn_roir]';
AllROI_m_sca_order=zeros(500,1);
AllROI_m_sca_order(All_label)=AllROI_m_sca;
Sig_sca=AllROI_m_sca_order(indnodes);
zerosind=find(Sig_sca==0);

Mean_rapid_slopecorr=Mean_rapid_slope';
Mean_rapid_slopecorr(zerosind)=[];
Sig_scacorr=Sig_sca;
Sig_scacorr(zerosind)=[];
%cov(zerosind,:)=[];
corrstats = corrstat(Mean_rapid_slopecorr,Sig_scacorr,[]);
corrplot(Mean_rapid_slopecorr,corrstats,Xlabels,Ylablels,[],[],[],[],[],[namef,'.tif'])