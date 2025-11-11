% using subset1 as an example
% For multicenter fitting, using multi_center_multimodelnlm_global.m instead of multimodelnlm_global.m
% For multicenter fitting, network atrriobutes should be harmonized (Adj_site_effect.m) before using multi_center_multimodelnlm_global.m
% 
statdir=[Pathnew,'Forgithub\Code'];
savedir=[Pathnew,'Forgithub'];
Pathnewplot=[Pathnew,'Forgithub'];
cd(Pathnewplot)

load([Pathnew,'Forgithub\Data\Gretna_results\Subset1\gretna_results_subset1.mat']);
load([Pathnew,'Forgithub\Data\Network\Subset1\Matrixallpar500FN_subset1.mat']);

Matrixallpar500FN=Matrixallpar500FN_subset1;
Age_gender=Age_gender_Brainsize_subset1(:,1:2);
Brainsz=Age_gender_Brainsize_subset1(:,3);
Eg=Eg_subset1;
Eloc=Eloc_subset1;
Ne=Ne_subset1;

[roi_num]=size(Matrixallpar500FN{1,1},2);
indall_select=[1:length(Brainsz)]';
plotornot=1;

cd(Pathnewplot)

%% Fig2 
Sv1=9.55;
Nname='NgE';
netname='par500FN';
Nodalmetrics=Ne;
opthub=[Pathnew,'Forgithub\Code\opt\nodes_opt.mat'];
nodal_stats(Nodalmetrics,Age_gender,Brainsz,Pathnewxsl,opthub,Pathnewplot,Sv1,Nname,netname,indall_select,1,1,1);
%nodal_stats2(Nodalmetrics,Age_gender,Brainsz,Pathnewxsl,opthub,Pathnewplot,Sv1,Nname,netname,indall_select,[Pathnewplot,'\indminNgEStrpar500FN0606.mat']);
opthub=[Pathnew,'Forgithub\Code\opt\nodes_slope_opt.mat'];
nodal_stats2(Pathnewxsl,opthub,Pathnewplot,Sv1,Nname,1,[Pathnewplot,'\indminNgEStr',netname,'.mat']);

