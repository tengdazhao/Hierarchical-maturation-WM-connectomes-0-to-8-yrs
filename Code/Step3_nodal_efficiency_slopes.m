load([Pathnew,'Forgithub\Data\Gretna_results\par500FN\NodalEfficiency\NodalEfficiency.mat']);
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

