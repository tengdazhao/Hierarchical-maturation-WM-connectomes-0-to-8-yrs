function [Ne_adj]=Adj_site_effect(NodalE_All,indall_select_All,Site_info,Brainsz_All,Age_gender_sp_All)
Ne_adj=NodalE_All;
NodalE_Alln=NodalE_All(indall_select_All,:);
Site_infon=Site_info(indall_select_All);
Brainsz_Alln=Brainsz_All(indall_select_All);
Age_Alln=Age_gender_sp_All(indall_select_All,1);
Gender_Alln=Age_gender_sp_All(indall_select_All,2);
Sp_Alln=Age_gender_sp_All(indall_select_All,3);
[nSubj, nRegion] = size(NodalE_Alln);
sites = categorical(Site_infon);
y_long = NodalE_Alln(:);
subID = repmat((1:nSubj)', nRegion,1);
regionID = repmat((1:nRegion), nSubj, 1);  % nSubj × nRegion
regionID = regionID(:);   
brainsize_long = repmat(Brainsz_Alln(:,1), nRegion,1);
Age_long = repmat(Age_Alln(:,1), nRegion,1);
Gender_long = repmat(Gender_Alln(:,1), nRegion,1);
Sp_long = repmat(Sp_Alln(:,1), nRegion,1);
sites_long = repmat(sites, nRegion,1);
sites_long= categorical(sites_long); 
brainsize_long = double(brainsize_long); 
y_long = double(y_long);  

T_long = table(y_long, brainsize_long, sites_long, subID, regionID,Age_long,Gender_long,Sp_long);

% remove the site effect
fit_global = fitlme(T_long, 'y_long ~ Sp_long + brainsize_long  +Gender_long  +(Age_long|sites_long)');


site_re = randomEffects(fit_global);
siteTbl = categories(sites);  
siteBias = zeros(nSubj,1);

for s = 1:numel(siteTbl)
    idx = sites == siteTbl{s};
    siteBias(idx) = site_re(s);  
end
Ne_adj(indall_select_All,:) = NodalE_All(indall_select_All,:) - siteBias;
