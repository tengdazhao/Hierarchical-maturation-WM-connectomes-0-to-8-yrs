SphereSurf=cell(2,1);
SphereSurf{1}='Q1-Q6_RelatedParcellation210.L.sphere.32k_fs_LR.surf.gii';
SphereSurf{2}='Q1-Q6_RelatedParcellation210.R.sphere.32k_fs_LR.surf.gii';
lhgii=['/HeLabData2/ztd/0~8/Myatlas/Par_500_fslr-lh.gii'];
rhgii=['/HeLabData2/ztd/0~8/Myatlas/Par_500_fslr-rh.gii'];
%The above dir should be replaced by customized dirs of surface and
%labeling fies

thissurfl = gifti(lhgii)
labelssurfl= thissurfl.cdata;
labelssurflqni = unique(labelssurfl(find(labelssurfl>0)));
thissurfr = gifti(rhgii)
labelssurfr= thissurfr.cdata;
labelssurfrqni = unique(labelssurfr(find(labelssurfr>0)));


labelssurfqni=[labelssurflqni; labelssurfrqni];
uqn_roi=unique(labelssurfqni);
NumNode=length(uqn_roi);


LV=gifti(lhgii); % FS_LR 32k
RV=gifti(rhgii);
NumHemi=size(LV.cdata, 1);

%% Spin Test
AAL500LR=[LV.cdata; RV.cdata];
%LabelLR=zeros(size(AAL500LR));

LabelLR=AAL500LR;

Label=cell(2,1);
Label{1}=LabelLR(1:NumHemi, 1);
Label{2}=LabelLR(NumHemi+1:end, 1);

NumPerm=10000; 
Perm_500_Results = zeros(NumNode,NumPerm);

parpool
parfor i=1:NumPerm
    OutFile=['spintests/results/Rand_' sprintf('%.5d', i)];
    fprintf('Performing %s\n', OutFile);
    GetRotateLabel(SphereSurf, Label, OutFile);
    
end

%     
% end
% 
% RROI_mean_scalar=[];
% uqn_roir=unique(labelssurfr);
% n_roi=length(uqn_roir)
% for i=1:n_roi
%     rightsca_n=rightscalar(find(labelssurfr==uqn_roir(i)));
%     RROI_mean_scalar(i)=nanmean(rightsca_n(:));
%    % RROI_mean_scalar(i)=mean(rightsca_n);
% end
% 
% AllROI_m_sca = [LROI_mean_scalar RROI_mean_scalar];
% All_label=[uqn_roil; uqn_roir]';

