rgs2 = dlmread('200mM_rgs.dat');

time = [7.11E-01
6.06E+00
1.14E+01
1.67E+01
2.21E+01
2.74E+01
3.28E+01
3.81E+01
4.35E+01
4.88E+01
5.42E+01
5.95E+01
6.49E+01
7.02E+01
7.55E+01
8.09E+01];

figure
hold on
plot(time,1.075*max(rgs2'),'-','color',[0.2 0.2 0.2])
plot(time,0.925*max(rgs2'),'-','color',[0.2 0.2 0.2])
patch([time', fliplr(time')],[1.075*max(rgs2'), fliplr(0.925*max(rgs2'))],[0.5 0.5 0.5])
plot(time,1.075*min(rgs2'),'-','color',[0.2 0.2 0.2])
plot(time,0.925*min(rgs2'),'-','color',[0.2 0.2 0.2])
patch([time', fliplr(time')],[1.075*min(rgs2'), fliplr(0.925*min(rgs2'))],[0.5 0.5 0.5])
alpha(0.3)
for i = 1:16
    plot(time(i)*ones(1,5),rgs2(i,:),'ok','markersize',4,'markerfacecolor','k')
    plot(time(i),mean(rgs2(i,:)),'ok','markersize',8,'markerfacecolor','r')
end


box on
format_figure
xlabel('time (ms)')
ylabel('R_G (Å)')

hl = hline(26.7)
xlim([0 max(time)+1])