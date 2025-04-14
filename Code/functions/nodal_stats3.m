function nodal_stats3(Pathnewxsl,opthub,Pathnewplot,Scalev1,Nname,pthretype,restulsmat);
result_name=Nname;
load(restulsmat);
FourModelsigthre=arrange_nodal_Slope2(indminNodalStrmdall,indminNodalStr,notenoughNodalStr,Pathnewxsl,opthub,Scalev1,Pathnewplot,pthretype,result_name);
end