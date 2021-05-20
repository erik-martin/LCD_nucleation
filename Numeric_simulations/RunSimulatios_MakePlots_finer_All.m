clear all
format compact

NA=[100,150,200,300,400,500,800,1000];
% Alldatanumber=17:-1:1;
Alldatanumber=501:-25:101;
load(sprintf('finedata_RU_%i_Rep_%i',101,1));
for i1=1:length(Alldatanumber)
    ML{i1}=sprintf('%i mM NaCL',NAx(Alldatanumber(i1)));
end
ML';
tmin=t_save(end);
tN=length(t_save);

for i3=1:length(Alldatanumber)
for i4=1:RepN
    datanumber=Alldatanumber(i3);
    load(sprintf('finedata_RU_%i_Rep_%i',datanumber,i4),'t_save')
    if length(t_save)>tN
        tN=length(t_save);
    end
end
end
AllTime=nan(length(Alldatanumber),tN,RepN);
% AllOmass=
% AllON=
% AllClust
% All_n_infty


% RepN=1;
for i3=1:length(Alldatanumber)
for i4=1:RepN
    
    % mass of clusters
    datanumber=Alldatanumber(i3);
    load(sprintf('finedata_RU_%i_Rep_%i',datanumber,i4))
% RepN=1;
    AllSmall(1:size(P_save,1),:,i3,i4)=P_save;
    n_infty=n_infty_0*sum(P_save(:,1:20).*repmat(1:20,size(P_save,1),1),2)/Ntot;
    AllTime(i3,1:length(t_save),i4)=t_save;
    if t_save(end)<tmin
        tmin=t_save(end);
    end
    Omass=nan(1,length(t_save));
    ON=nan(1,length(t_save));
    ClusterThreshold=20;
    for t=1:size(O_save,1)
        Omass(t)=sum(O_save(t,O_save(t,:)>ClusterThreshold));
        ON(t)   =sum(O_save(t,:)>ClusterThreshold);
    end
    AllOmass(i3,1:length(t_save),i4)=Omass;
    AllON(i3,1:length(t_save),i4)=ON;
    AllClust(1:length(t_save),1:size(O_save,2),i3,i4)=O_save;
    
    HalfPoint(i3,i4)=t_save(find(Omass>Omass(end)/2,1));
    
    if i3==10 && i4==3
        HalfPoint(i3,i4)=t_save(find(Omass>Omass(end-5)/2,1));
    end
    
    All_n_infty(i3,1:length(t_save),i4)=n_infty;
    All_c_p(i3)=c_p;
%     [n_m,n_p]
    All_n_p(i3)=n_p;
    
end
end
AllNA=NAx(Alldatanumber);
MaxAllTime=max(max(max(AllTime)));
MaxAllTime=.035;



%% Plot mass of clusters
% tr=1:size(AllTime,2);
% figure(1)
% mc=winter(size(AllTime,1));
% % for i2=1:size(AllTime,3)
% for i1=1:size(AllTime,1)
%     plot(mean(AllTime(i1,tr,:),3),mean(AllOmass(i1,tr,:),3),'linewidth',3,'color',mc(i1,:))
%     hold on
% end
% % end
% hold off
% xlabel('time')
% xlim([0,tmin])
% ylabel(['Mass of clusters >',num2str(ClusterThreshold)])
% set(gca,'fontsize',16)
% legend(ML,'location','eastoutside')
% saveas(gca,['C_AllClusterMass_Mean'],'epsc')

% tr=1:size(AllTime,2);
% figure(1)
% plot(mean(AllTime(:,tr,:),3)',mean(AllOmass(:,tr,:),3)','linewidth',3)
% xlabel('time')
% xlim([0,tmin])
% ylabel(['Mass of clusters >',num2str(ClusterThreshold)])
% set(gca,'fontsize',16)
% % legend(ML,'location','southeast')
% saveas(gca,['C_AllClusterMass_Mean'],'epsc')

% tr=1:size(AllTime,2);
% figure(2)
% mc=winter(size(AllTime,1));
% for i2=1:size(AllTime,3)
% for i1=1:size(AllTime,1)
%     plot(squeeze(AllTime(i1,tr,i2)),squeeze(AllOmass(i1,tr,i2)),'linewidth',3,'color',mc(i1,:))
%     hold on
% end
% end
% hold off
% xlabel('time')
% xlim([0,tmin])
% ylabel(['Mass of clusters >',num2str(ClusterThreshold)])
% set(gca,'fontsize',16)
% legend(ML,'location','eastoutside')
% saveas(gca,['C_AllClusterMass_Reps'],'epsc')
% return

%% plot dilute concentration
% figure(3)
% plot(AllTime(:,tr)',All_n_infty'/RescaleC,'linewidth',3);
% set(gca,'ColorOrderIndex',1)
% hold on
% plot([0,t_save(end)],[All_c_p;All_c_p],'--','linewidth',3)
% hold off
% xlim([0,tmax])
% legend(ML,'location','northeast')



% return

%% plot halfway point (regular droplet growth)
hp=mean(HalfPoint(end-2,:));
ss=400-cSat_Allx(Alldatanumber(end-2));
rHalfPoint=HalfPoint(1:end-2,:);
SuperSat=400-cSat_Allx(Alldatanumber(1:end-2));
% pf=polyfit(SuperSat,mean(log(rHalfPoint-B),2),1);
% a0=[pf(1),pf(2),.003];
% F=@(a,xdata)exp(a(1)*xdata+a(2)) +a(3);

a0=[.0031*300^(1/5)*exp(-(1000/300)^2),1000];
F=@(a,xdata)a(1)*xdata.^(-1/5).*exp((a(2)./xdata).^2);
logF=@(a,xdata)log(F(a,xdata));
x=SuperSat;
y=mean(rHalfPoint,2)';
opts = optimset('Display','off');
[a,resnorm,~,exitflag,output] = lsqcurvefit(logF,a0,x,log(y),[ ],[ ],opts);
figure(4)
xf=360:-.01:400-cSat_Allx(Alldatanumber(end));
xf=360:-.01:200;
h0=plot(xf,F(a,xf),'linewidth',3);
hold on
h1=plot(x(1),rHalfPoint(1),'ok','markersize',6);
h1_all=plot(x,rHalfPoint,'ok','markersize',6)
% plot(ss,hp,'.r','markersize',16)
SuperSat_notfinished=400-cSat_Allx(Alldatanumber(end-1:end));
h2=plot(SuperSat_notfinished,F(a,SuperSat_notfinished),'kd','markersize',6);
hold off
ylabel('Time to half completion (seconds)')
xlabel('Supersaturation: C-C_{sat} (\muM)')
set(h1_all, 'markerfacecolor', get(h1, 'color'));
set(h2, 'markerfacecolor', get(h2, 'color'));
set(h1, 'markerfacecolor', get(h1, 'color'));
legend([h1,h0,h2],'individual simulation','fit','extrapolated time')
set(gca,'fontsize',16)
xlim([230,345])
saveas(gca,'C_LagTime_regular_Lin','epsc')
set(gca,'yscale','log')
saveas(gca,'C_LagTime_regular_Log','epsc')



%% plot halfway point
hp=mean(HalfPoint(end-2,:));
ss=400-cSat_Allx(Alldatanumber(end-2));
rHalfPoint=HalfPoint(1:end-2,:);
SuperSat=400-cSat_Allx(Alldatanumber(1:end-2));
% pf=polyfit(SuperSat,mean(log(rHalfPoint-B),2),1);
% a0=[pf(1),pf(2),.003];
% F=@(a,xdata)exp(a(1)*xdata+a(2)) +a(3);

a0=[.0031*exp(-(1000/300)^2),1000];
F=@(a,xdata)a(1).*exp((a(2)./xdata).^2);
logF=@(a,xdata)log(F(a,xdata));
x=SuperSat;
y=mean(rHalfPoint,2)';
opts = optimset('Display','off');
[a,resnorm,~,exitflag,output] = lsqcurvefit(logF,a0,x,log(y),[ ],[ ],opts);
figure(5)
xf=360:-.01:400-cSat_Allx(Alldatanumber(end));
xf=360:-.01:200;
h0=plot(xf,F(a,xf),'linewidth',3);
hold on
h1=plot(x(1),rHalfPoint(1),'ok','markersize',6);
h1_all=plot(x,rHalfPoint,'ok','markersize',6)
% plot(ss,hp,'.r','markersize',16)
SuperSat_notfinished=400-cSat_Allx(Alldatanumber(end-1:end));
h2=plot(SuperSat_notfinished,F(a,SuperSat_notfinished),'kd','markersize',6);
hold off
ylabel('Time to half completion (seconds)')
xlabel('Supersaturation: C-C_{sat} (\muM)')
set(h1_all, 'markerfacecolor', get(h1, 'color'));
set(h2, 'markerfacecolor', get(h2, 'color'));
set(h1, 'markerfacecolor', get(h1, 'color'));
legend([h1,h0,h2],'individual simulation','fit','extrapolated time')
set(gca,'fontsize',16)
xlim([230,345])
saveas(gca,'C_LagTime_stagnant_Lin','epsc')
set(gca,'yscale','log')
saveas(gca,'C_LagTime_stagnant_Log','epsc')
return
%% Plot All Clusters
% figure(5)
% for i1=1:size(AllClust,3)
% for i2=1:size(AllClust,4)
%     t=AllTime(i1,:,i2);
%     y=sort(AllClust(:,:,i1,i2)','descend')';
%     ff=sum(y>0,1)>0;
%     
%     plot(AllTime(i1,:,i2),y(:,ff))
% %     plot(AllTime(i1,:,i2),max(AllClust(:,:,i1,i2)')')
%     xlabel('time (seconds)')
%     ylabel('Cluster size (number of proteins)')
%     set(gca,'fontsize',16)
%     
%     
%     ym=max(max(AllClust(:,:,i1,i2)));
%     yl=ceil(ym/10^floor(log10(ym)))*10^floor(log10(ym));
%     ylim([0,yl])
%     text(.002,yl*.9,sprintf('%i mM NaCL',AllNA(i1)),'fontsize',16)
%     print('-depsc2','-painters',sprintf('E_Clusters_C_%i_R_%i',i1,i2))
% end
% end
% 
% return

% for i=1:20
%     leg(i)={sprintf('cluster of size %i',i)};
% end
% figure(1)
% for i1=1:size(AllClust,3)
% for i2=1:size(AllClust,4)
%     plot(AllTime(i1,:,i2),AllSmall(:,1:20,i1,i2))
%     legend(leg,'location','eastoutside')
%     set(gca,'yscale','log')
%     set(gca,'fontsize',16)
%     xlabel('time (seconds)')
%     ylabel('Number of clusters')
%     print('-depsc2','-painters',sprintf('E_SmallClusters_C_%i_R_%i',i1,i2))
% end
% end

% tl=[10.^(0:6)];
% v=num2str(tl');
% figure(1)
% n=0;
% AllTime(:,end+1:end+2,:)=MaxAllTime;
% for i1=1:size(AllClust,3)
% for i2=1:size(AllClust,4)
%     y=[AllSmall(:,1:20,i1,i2),zeros(size(AllSmall,1),1)];
%     y=[y;zeros(2,size(y,2))];
%     h=pcolor(AllTime(i1,:,i2),1:21,log10(y)');
%     set(h, 'EdgeColor', 'none');
%     set(gca,'fontsize',16)
%     xlabel('time (seconds)')
%     ylabel('Cluster size')
%     
%     h2=colorbar;
%     cmap=colormap;
%     cmap(1,:)=1;
%     colormap(cmap);
%     
%     caxis([-7/length(cmap),6])
%     set(h2,'TickLabels',v);
%     cmap=colormap;
%     cmap(1,:)=1;
%     colormap(cmap);
%     
%     print('-depsc2','-painters',sprintf('F_HMSmallClusters_C_%i_R_%i',i1,i2))
% end
% end
% return

%% Plot All Clusters
figure(5)
for i1=10%1:size(AllClust,3)
for i2=3%1:size(AllClust,4)
    t=AllTime(i1,:,i2);
    y=sort(AllClust(:,:,i1,i2)','descend')';
    ff=sum(y>0,1)>0;
    figure(5)
    if i1==10 && i2==3
        plot(t(1:997),y(1:997,ff))
    else
        plot(t,y(:,ff))
    end
    xlabel('time (seconds)')
    ylabel('Cluster size (number of proteins)')
    set(gca,'fontsize',16)
    xlim([0,MaxAllTime])
    
    ym=max(max(AllClust(:,:,i1,i2)));
    yl=ceil(ym/10^floor(log10(ym)))*10^floor(log10(ym));
    ylim([-.5,yl])
    text(.002,yl*.9,sprintf('%i mM NaCL',AllNA(i1)),'fontsize',16)
%     saveas(gca,sprintf('D_Clusters_C_%i_R_%i',i1,i2),'epsc')
    print('-depsc2','-painters',sprintf('F_Clusters_C_%i_R_%i',i1,i2))
end
end

return
figure(2)
for i1=1:size(AllClust,3)
    for i2=1:size(AllClust,4)
        plot(sum(AllSmall(:,:,i1,i2),1)/sum(~isnan(AllTime(i1,:,i2))),'.','markersize',16)
        hold on
        set(gca,'yscale','log')
    end
    hold off
    xlabel('Cluster size')
    ylabel('Average number of clusters observed')
    legend('Run 1','Run 2','Run 3')
    set(gca,'fontsize',16)
    print('-depsc2','-painters',sprintf('E_ClusterDensity_C_%i_R_%i',i1,i2))
end



return

%% Nucleation size
All_n_p=repmat(All_n_p',1,size(All_n_infty,2),size(All_n_infty,3));
NN= 4/3*pi*n_m*(2*gamma*All_n_p./(n_m*(All_n_infty-All_n_p))).^3;
figure(5)
tr=1:size(AllTime,2);
mc=winter(size(AllTime,1));
for i1=1:size(AllTime,1)
    plot(mean(AllTime(i1,tr,:),3),mean(NN(i1,tr,:),3),'linewidth',3,'color',mc(i1,:))
    hold on
end
hold off
xlabel('time (seconds)')
xlim([0,tmin])
ylim([0,1000])
ylabel('Nucleation size (number of proteins)')
set(gca,'fontsize',16)
% legend(ML,'location','eastoutside')
saveas(gca,'C_NucleationSize_mean','epsc')


% All_n_p=repmat(All_n_p',1,size(All_n_infty,2),size(All_n_infty,3));
NN= 4/3*pi*n_m*(2*gamma*All_n_p./(n_m*(All_n_infty-All_n_p))).^3;
figure(5)
tr=1:size(AllTime,2);
mc=winter(size(AllTime,1));
for i2=1:size(AllTime,3)
for i1=1:size(AllTime,1)
    plot(squeeze(AllTime(i1,tr,i2)),squeeze(NN(i1,tr,i2)),'linewidth',3,'color',mc(i1,:))
    hold on
end
end
hold off
xlabel('time (seconds)')
xlim([0,tmin])
ylim([0,1000])
ylabel('Nucleation size (number of proteins)')
set(gca,'fontsize',16)
legend(ML,'location','eastoutside')
saveas(gca,'C_NucleationSize_Reps','epsc')