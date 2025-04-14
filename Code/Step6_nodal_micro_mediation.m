%% Fig5 
Sv1=9.55;
Nname='NgE';
netname='indminNgEStrpar500FN';
Nodalmetrics=Ne;
opthub=[Pathnew,'0-8years\stats\FinalFigs\nodes_slope_opt.mat'];
nodal_stats3(Pathnewxsl,opthub,Pathnewplot,Sv1,Nname,1,[Pathnewplot,'\',netname,'.mat']);

