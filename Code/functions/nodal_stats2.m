function nodal_stats2(Nodalmetric,Pathnewxsl,opthub,Pathnewplot,Scalev1,Nname,pthretype,restulsmat);
result_name=Nname;
load(restulsmat);
FourModelsigthre=arrange_nodal_Slope(Nodalmetric,indminNodalStrmdall,indminNodalStr,notenoughNodalStr,Pathnewxsl,opthub,Scalev1,Pathnewplot,pthretype,result_name);
end