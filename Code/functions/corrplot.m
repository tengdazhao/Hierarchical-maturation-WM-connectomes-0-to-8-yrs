function corrplot(data2,corr_stats,xlabeln,ylabeln,xtk,xtklabel,ytk,ytklabel,ppp,fname)
ax = min(data2):(max(data2)-min(data2))/100:max(data2);
p = corr_stats{1,4};
r = corr_stats{1,1};

figure
N=0;
pwhile = p;
if(p>0)
while(pwhile<1)
    pwhile = pwhile*10;
    N=N+1;
end

pstr = num2str(pwhile);
pstr = pstr(1:3);

pstr = ['p = ',pstr, ' x 10^{\bf\fontsize{28}-' num2str(N) '}'];
% pstr = ['p < 10^{\bf\fontsize{32}-20}'];
else    
pstr = ['p = 0;'];
end
if ~isempty(ppp)
    if(ppp)>0.001
    pstr = ['p = ' num2str(ppp)];
    elseif(ppp)<=0.001
        pstr = ['p < ' num2str(ppp)];
    end
end
rkeep = vpa(r,2);
% if(mod(r*100,10)==0)
if(length(char(rkeep))==3)
 rstr = ['r = ',char(rkeep),'0;'];
 else
rstr = ['r = ',char(rkeep),';'];
 end
sca = corr_stats{1,3};

whichstats={'rsquare','tstat','mse','r'};
%linear
stats_linear = regstats(sca,data2,'linear',whichstats);
beta = stats_linear.tstat.beta;
% data_fitted = beta(1)+beta(2)*data2+stats_linear.r;
line = (beta(1)+ beta(2)*ax);
% figure;
axes('LineWidth',5,'FontWeight','bold','FontSize',30);
axis([min(data2)-(max(data2)-min(data2))/10 max(data2)+(max(data2)-min(data2))/12 min(sca)-(max(sca)-min(sca))/6 max(sca)+(max(sca)-min(sca))/4])
%axis([min(data2)-(max(data2)-min(data2))/10 max(data2)+(max(data2)-min(data2))/12 -4 55]);
set(gcf,'Color',[1 1 1]) ;
set(gcf,'position', [1 30 950 715]);
set(gca,'position',[ 0.23    0.22   0.71    0.74]);
hold on;
plot(data2,sca,'o','LineWidth',2,'MarkerFaceColor',[0 0 0],'MarkerSize',12,'MarkerEdgeColor',[0 0 0])
plot(ax,line,'k','LineWidth',8,'LineStyle','-','color',[216/255    41/255    0/255])
xlabel(xlabeln);
ylabel(ylabeln);
% set(gca,'ytick',[str2double(char(vpa(max(sca)-(max(sca)-round(mean(max(sca)+min(sca))/2))*2,2))),str2double(char(vpa(round((max(sca)-min(sca))/2+min(sca)),2))),str2double(char(vpa(max(sca),2)))]); 
% set(gca,'xtick',[str2double(char(vpa(max(data2)-(max(data2)-round(mean(max(data2)+min(data2))/2))*2,2))),str2double(char(vpa(round((max(data2)-min(data2))/2+min(data2)),2))),str2double(char(vpa(max(data2),2)))]); 
tempstry = ((max(sca)-min(sca))/2+min(sca));
fwei = floor(log10(tempstry));
if(fwei<2);
    fwei = abs(fwei);   
tempstry = tempstry*(fwei+1)*10;
round_tempstry = round(tempstry);
round_tempstry = round_tempstry/((fwei+1)*10);
else
round_tempstry = round(tempstry);
end
    
tempstrx = ((max(data2)-min(data2))/2+min(data2));
fwei = floor(log10(tempstrx));
if(fwei<2);
    fwei = abs(fwei);   
tempstrx = tempstrx*(fwei+1)*10;
round_tempstrx = round(tempstrx);
round_tempstrx = round_tempstrx/((fwei+1)*10);
else
round_tempstrx = round(tempstrx);
end
    
% set(gca,'ytick',[str2double(char(vpa(min(sca),2))),str2double(char(vpa(round((max(sca)-min(sca))/2+min(sca)),2))),str2double(char(vpa(max(sca),2)))]); 
% set(gca,'xtick',[str2double(char(vpa(min(data2),2))),str2double(char(vpa(round((max(data2)-min(data2))/2+min(data2)),2))),str2double(char(vpa(max(data2),2)))]); 

set(gca,'ytick',[str2double(char(vpa(min(sca),3))),str2double(char(vpa(round_tempstry,3))),str2double(char(vpa(max(sca),3)))]); 
set(gca,'xtick',[str2double(char(vpa(min(data2),3))),str2double(char(vpa(round_tempstrx,3))),str2double(char(vpa(max(data2),3)))]); 
% 
if ~isempty(ytk)
set(gca,'ytick',ytk,'yticklabel',ytklabel); 
set(gca,'xtick',xtk,'xticklabel',xtklabel); 
end
hold on;
% set(gca,'ytick',[20 40 60 80]);
text(min(data2)+(max(data2)-min(data2))/2.6,max(sca)+(max(sca)-min(sca))/5, pstr,'FontWeight','bold','FontSize',40);
text(min(data2)-(max(data2)-min(data2))/22,max(sca)+(max(sca)-min(sca))/5.5, rstr,'FontWeight','bold','FontSize',40);
hold on;
set(gca,'fontsize',40,'LineWidth',5,'FontWeight','bold');
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf,'-dtiff','-r500',[fname '.tif']);  