function [indmin mdall notenough]=multimodelnlm_global(Age_gender,Brainsz,data_r,ylabelname,Plotornot,indall_select,betarealinitial)
for i=1:size(indall_select,2);
Age_gender_new=Age_gender(indall_select(:,i),:);
Brainsz_new=Brainsz(indall_select(:,i));
data_r_new=data_r(indall_select(:,i));
[indmin(:,:,i) mdall{i} notenough{i}]=multimodelnlm_global_novalid(Age_gender_new,Brainsz_new,data_r_new,['Split',num2str(i),'  ',ylabelname],Plotornot,betarealinitial);
end
end
