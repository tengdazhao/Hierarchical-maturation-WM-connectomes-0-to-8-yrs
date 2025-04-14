function nodal_stats(NodalStr,Age_gender,Brainsz,Pathnewxsl,opthub,Pathnewplot,Scalev1,Nname,netname,indall_select,betainitial,plotornot,fdrbon);
% %sort
Strsort=[];indstr=[];indstrre=[];re=[];
roi_num=size(NodalStr,2);

indminNodalStr=[];
indminNodalStrmdall=[];
for i=1:roi_num%60
    Nodal=NodalStr(:,i);
    Namen=[Nname, num2str(i)];
myflag = true;
trynum=0;
while myflag
 trynum=trynum+1;
    try
        if(trynum>100)
   [indminNodalStr{i}.ModelselectRestults indminNodalStrmdall{i} notenoughNodalStr{i}]=multimodelnlm_globalonlylinear(Age_gender,Brainsz,Nodal,Namen,0,indall_select,betainitial);
    myflag = false;
    break
         end
        
       if(trynum<=20)
%
        [indminNodalStr{i}.ModelselectRestults indminNodalStrmdall{i} notenoughNodalStr{i}]=multimodelnlm_global(Age_gender,Brainsz,Nodal,Namen,0,indall_select,betainitial);
      
       else(trynum>20&&trynum<100)
        [indminNodalStr{i}.ModelselectRestults indminNodalStrmdall{i} notenoughNodalStr{i}]=multimodelnlm_global(Age_gender,Brainsz,Nodal,Namen,0,indall_select,0);
       end
    myflag = false;
    
    end
end
i
end


fileforsavename=[Pathnewplot,'indmin',Nname,'Str',netname,'.mat'];
save(fileforsavename,'indminNodalStr','indminNodalStrmdall','notenoughNodalStr');

result_name=[Nname,'Str'];
FourModelsigthre=arrange_nodal(indminNodalStrmdall,indminNodalStr,notenoughNodalStr,Pathnewxsl,opthub,Scalev1,Pathnewplot,result_name,plotornot,fdrbon);

fileforsavename=[Pathnewplot,'FourModelsigthre',Nname,'Str',netname,'.mat'];
save(fileforsavename,'FourModelsigthre');

end