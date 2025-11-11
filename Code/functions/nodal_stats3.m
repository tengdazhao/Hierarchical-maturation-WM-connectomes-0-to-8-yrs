function nodal_stats3(Nodal_metrics,Pathnewxsl,opthub,Pathnewplot,Scalev1,Nname,pthretype,restulsmat);
result_name=Nname;
load(restulsmat);
FourModelsigthre=arrange_nodal_Slope2(Nodal_metrics,indminNodalStrmdall,indminNodalStr,notenoughNodalStr,Pathnewxsl,opthub,Scalev1,Pathnewplot,pthretype,result_name);
end