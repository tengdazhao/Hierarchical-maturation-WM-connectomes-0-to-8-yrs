function [indmin mdall notenough]=multicenter_multimodelnlm_global(Age_gender,Brainsz,data_r,ylabelname,Plotornot,indall_select,betarealinitial,Site_info)
for i=1:size(indall_select,2);
Age_gender_new=Age_gender(indall_select(:,i),:);
Brainsz_new=Brainsz(indall_select(:,i));
data_r_new=data_r(indall_select(:,i));
Site_info_new=Site_info(indall_select(:,i));
%Site_info_new=Site_info;
[indmin(:,:,i) mdall{i} notenough{i}]=multicenter_multimodelnlm_global_novalid(Age_gender_new,Brainsz_new,data_r_new,['Split',num2str(i),'  ',ylabelname],Plotornot,betarealinitial,Site_info_new);
end
end
