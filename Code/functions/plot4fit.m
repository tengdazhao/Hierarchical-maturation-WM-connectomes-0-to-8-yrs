function plot4fit(Age,Data_rremain,fitresult,name,ylabelname,pvalue,rsq)
N=0;
pwhile = pvalue;
if(pvalue>0)
while(pwhile<1)
    pwhile = pwhile*10;
    N=N+1;
end

pstr = num2str(pwhile);
if(length(pstr)<3)
pstr = pstr;
else
pstr = pstr(1:3);
end
pstr = ['p = ',pstr, ' x 10^{\bf\fontsize{28}-' num2str(N) '}'];
% pstr = ['p < 10^{\bf\fontsize{32}-20}'];
else    
pstr = ['p = 0;'];
end
ppp=[];% manual set p value
if ~isempty(ppp)
    if(ppp)>0.001
    pstr = ['p = ' num2str(ppp)];
    elseif(ppp)<=0.001
        pstr = ['p < ' num2str(ppp)];
    end
end
rkeep = vpa(rsq,2);
% if(mod(r*100,10)==0)
if(length(char(rkeep))==3)
 rstr = ['r = ',char(rkeep),'0;'];
 else
rstr = ['r = ',char(rkeep),';'];
 end


close all
figure;
yData=Data_rremain;
xData=Age;
axes('LineWidth',5,'FontWeight','bold','FontSize',48);
axis([min(xData)-(max(xData)-min(xData))/10 max(xData)+(max(xData)-min(xData))/12 min(yData)-(max(yData)-min(yData))/5 max(yData)+(max(yData)-min(yData))/5]);
hold on;

%h=plot(fitresult,Age,Data_rremain, 'predobs');
h=plot(fitresult,Age,Data_rremain, 'predfun');

%set(h(1),'LineWidth',1.2,'Marker','o','MarkerFaceColor',[42/255    158/255    237/255],'MarkerSize',12,'MarkerEdgeColor',[0 0 0]);
%set(h(1),'LineWidth',0.6,'Marker','o','MarkerFaceColor',[128/255    198/255    255/255],'MarkerSize',9,'MarkerEdgeColor',[0 0 0]);
%set(h(1),'LineWidth',1.2,'Marker','o','MarkerFaceColor',[32/255    138/255    217/255],'MarkerSize',12,'MarkerEdgeColor',[0 0 0]);
set(h(1),'LineWidth',0.6,'Marker','o','MarkerFaceColor',[20/255    158/255    205/255],'MarkerSize',9,'MarkerEdgeColor',[0 0 0]);
set(h(2),'LineWidth',4.5,'LineStyle','-','color',[255/255 0 0]);
%set(h(3),'LineWidth',2,'LineStyle','-','color',[1    150/255    150/255]);
%set(h(4),'LineWidth',2,'LineStyle','-','color',[1    150/255    150/255]);
set(h(3),'LineWidth',0.1,'LineStyle','-','color',[255/255     224/255    200/255]);
set(h(4),'LineWidth',0.1,'LineStyle','-','color',[255/255     224/255    200/255]);
fill([h(3).XData,fliplr(h(3).XData)],[h(3).YData,fliplr(h(4).YData)],[255/255     169/255    160/255],'EdgeColor','none','FaceAlpha',0.5);%»­Ìî³äÑÕÉ«
%h=plot(fitresult,Age,Data_rremain);
plot(h(2).XData,h(2).YData,'LineWidth',4.5,'LineStyle','-','color',[255/255 0 0]);
% set(h(3),'Visible','off');
% set(h(4),'Visible','off');

% xp = min(peak(4:5)):0.01:max(peak(4:5));
% yp = zeros(1,length(xp))+(min(yData)-0.25*(max(yData)-min(yData)));
% % plot(xp,yp,'LineWidth',9,'LineStyle','-','color',[0.98   0.55    0]);
% % yp2 = (min(yData)-0.35*(max(yData)-min(yData))):0.02*0.2*(max(yData)-min(yData)):(min(yData)-0.25*(max(yData)-min(yData)));
% % plot(zeros(1,length(yp2))+peak(5),yp2,'LineWidth',6,'LineStyle','-','color',[0.98   0.55    0]);
% % plot(zeros(1,length(yp2))+peak(4),yp2,'LineWidth',6,'LineStyle','-','color',[0.98   0.55    0]);
% plot(xp,yp,'LineWidth',9,'LineStyle','-','color',[255/255   204/255    0]);
% yp2 = (min(yData)-0.3*(max(yData)-min(yData))):0.02*0.1*(max(yData)-min(yData)):(min(yData)-0.2*(max(yData)-min(yData)));
% plot(zeros(1,length(yp2))+peak(5),yp2,'LineWidth',6,'LineStyle','-','color',[255/255   204/255    0]);
% plot(zeros(1,length(yp2))+peak(4),yp2,'LineWidth',6,'LineStyle','-','color',[255/255   204/255    0]);
% Label axes
% set(gca,'ytick',[str2double(char(vpa(min(yData),2))),str2double(char(vpa(round((max(yData)-min(yData))/2+min(yData)),2))),str2double(char(vpa(max(yData),2)))]); 
posorneg=log10((max(yData)-(min(yData)))/2);
if(posorneg>=0)
Nz=10^(round(posorneg)-1);


elseif(posorneg<0)
Nz=10^(round(posorneg)); 


end
yData=yData/Nz;
Yt=[double(round(vpa(max(yData)-(max(yData)-round(mean(max(yData)+min(yData))/2))*2,2))*Nz),double(round(vpa(round((max(yData)-min(yData))/2+min(yData)),2))*Nz),double(round(vpa(max(yData),2))*Nz)];
yData=yData*Nz;
cn=1;
while((length(Yt)-length(unique(Yt)))) 
Nz=10^(round(posorneg)-cn); 
cn=cn+1;
yData=yData/Nz;
Yt=[double(round(vpa(max(yData)-(max(yData)-round(mean(max(yData)+min(yData))/2))*2,2))*Nz),double(round(vpa(round((max(yData)-min(yData))/2+min(yData)),2))*Nz),double(round(vpa(max(yData),2))*Nz)];
yData=yData*Nz;
end
set(gca,'ytick',Yt);

text(min(xData)+(max(xData)-min(xData))/2.6,max(yData)+(max(yData)-min(yData))/5, pstr,'FontWeight','bold','FontSize',40);
text(min(xData)-(max(xData)-min(xData))/22,max(yData)+(max(yData)-min(yData))/5.5, rstr,'FontWeight','bold','FontSize',40);

% set(gca,'XTick',[10 20 30 40 50 60 70 80]);
% set(gca,'XtickLabel',['  ';'20';'  ';'40';'  ';'60';'  ';'80']);
xlabel( 'Age','FontSize',52);
ylabel(ylabelname,'FontSize',52);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 70 970 800]);
set(gca,'position',[ 0.25    0.24   0.71  0.71]);
set(gcf, 'PaperPositionMode', 'auto');
legend('off');
print(gcf,'-dtiff','-r300',[name]);
close all
end
