%% Fig5 
% set working paths
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


Sv1=9.55;
Nname='NgE';
netname='indminNgEStrpar500FN';
Nodalmetrics=Ne;

opthub=[Pathnew,'Forgithub\Code\opt\nodes_slope_opt.mat'];
nodal_stats3(Pathnewxsl,opthub,Pathnewplot,Sv1,Nname,1,[Pathnewplot,'\',netname,'.mat']);

