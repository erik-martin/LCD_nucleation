% use this script to make plots of Rg versus Time for A1-LCD and Aro
% Depleted LCD
% the function "format_figure" contains a set of parameters to modify the
% figure. It can be commented out. 

% Rg/ I0 from MFF
% directory:
cd /A1LCD/mff
Rgs = dlmread('Rg.dat');
RgsE = dlmread('RgE.dat');
I0s = dlmread('I0.dat');
I0sE = dlmread('I0E.dat');

tm = [6.88E+01
1.46E+02
2.09E+02
3.02E+02
4.62E+02
6.46E+02
9.26E+02
1.21E+03
1.53E+03
1.90E+03
2.31E+03
2.77E+03
3.27E+03
3.81E+03
4.39E+03
5.02E+03
5.70E+03
6.42E+03
7.18E+03
8.25E+03
9.29E+03
1.03E+04
1.19E+04
1.33E+04
1.47E+04
1.61E+04
1.74E+04];


tm = tm/1000000


%plotting
sz = [[440   378   800   350]];
fh= figure('position',sz)

set(fh,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);

yyaxis left
hold on


Rgs2E = Rgs.^2 .* 2 .* (RgsE./Rgs);
for i = 1:length(tm);
    plot([tm(i) tm(i)], [Rgs(i).^2+Rgs2E(i) Rgs(i).^2-Rgs2E(i)],'-k','linewidth',0.6)
end
plot(tm,Rgs.^2,'ok','markersize',8,'markerfacecolor','k');
ylabel('R_G^2')



eqsum = @(c,d)c(1).*exp(c(2).*d) + c(3).*(1-exp(c(4).*d)) + c(5);

d0 = [5.4206 -545.7054 20.0000  -36.8653   19.3862];
%d0 = [286.8   -405.0022 5.1345e+04 0.4829 -5.0952e+04];
c0min = [0 -inf 0]; 
c0max = [20 inf inf];
fit3 = lsqcurvefit(eqsum,d0,tm,Rgs.^2);
ymodel3 = eqsum(fit3,xmodel);

plot(xmodel,ymodel3,'-k')

%ylim([21 30])
xlim([0 0.018])
yyaxis right
hold on


for i = 1:length(tm);
    plot([tm(i) tm(i)], [I0s(i)+I0sE(i) I0s(i)-I0sE(i)],'-r','linewidth',0.6)
end
plot(tm,I0s,'o','markersize',8,'markerfacecolor','r');
hline(0.2153)
ylabel('I_0')
ylim([0 1])
xlabel('time (sec)')
box on
format_figure



%% aromatic depleted data
% Rg/ I0 from MFF for aro2 data

cd AroDepletedLCD/mff
Rgs = dlmread('Rg.dat');
RgsE = dlmread('RgE.dat');
I0s = dlmread('I0.dat');
I0sE = dlmread('I0E.dat');

tm = [0.00096875
0.00184375
0.00271875
0.00376875
0.00464375
0.00551875
0.00639375
0.00726875
0.00814375
0.00901875
0.00989375
0.0107687
0.0116437
0.0125187
0.0133937
0.0142687
0.0151437
0.0160187
0.0168937
0.0177687
0.0186437
0.0195187
0.0200437];



% plotting

sz = [[440   378   800   350]];
fh= figure('position',sz)

set(fh,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);



yyaxis left
hold on


Rgs2E = Rgs.^2 .* 2 .* (RgsE./Rgs);
for i = 1:length(tm);
    plot([tm(i) tm(i)], [Rgs(i).^2+Rgs2E(i) Rgs(i).^2-Rgs2E(i)],'-k','linewidth',0.6)
end
plot(tm,Rgs.^2,'ok','markersize',8,'markerfacecolor','w');
ylabel('R_G^2')


expf = @(c,d) c(1).*exp(c(2).*d) + c(3);
expf2 = @(c,d) c(1).*(1-exp(c(2).*d)) + c(3);
sigf = @(a,b)a(1) - a(2)./(1+ exp(-a(3).*(b-a(4))));


a0 = [.001 -.1 30];
a0min = [1 -inf -inf]; 
a0max = [inf 0 inf];
fit1 = lsqcurvefit(expf,a0,tm(1:15),Rgs(1:15),a0min,a0max);
xmodel = 0:0.000001:0.018;
ymodel1 = expf(fit1,xmodel);
plot(xmodel,ymodel1.^2,'--k')

c0 = [1 .1 16];
c0min = [0 -inf 0]; 
c0max = [20 inf inf];
fit2 = lsqcurvefit(expf2,c0,tm(15:end),Rgs(15:end),c0min,c0max);
ymodel2 = expf2(fit2,xmodel);



xlim([0 0.018])
ylim([600 1100])
yyaxis right
hold on


for i = 1:length(tm);
    plot([tm(i) tm(i)], [I0s(i)+I0sE(i) I0s(i)-I0sE(i)],'-r','linewidth',0.6)
end
plot(tm,I0s,'or','markersize',8,'markerfacecolor','w');
hline(mean(I0s))
ylabel('I_0')
ylim([0.75 1.8])
xlabel('time (sec)')
box on



xlabel('time (sec)')
box on
format_figure


