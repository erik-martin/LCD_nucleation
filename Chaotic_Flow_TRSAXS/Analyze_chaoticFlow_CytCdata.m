%plot Cyt C chaotic flow Rg and I0 data

cd CytC
c2_1 = importdata('12_2019_avg_cytc2_avg_Rgs.csv',',');
c2_1 = c2_1.data;
c2_2 = importdata('08_2019_avg_cytc_avg_Rgs.csv',',');
c2_2 = c2_2.data;

figure
hold on

%prop error to square:
c2_1e = c2_1(:,3).^2 .* 2 .* (c2_1(:,4)./c2_1(:,3));
for i = 1:length(c2_1(:,2));
    plot([c2_1(i,2) c2_1(i,2)], [c2_1(i,3).^2+c2_1e(i) c2_1(i,3).^2-c2_1e(i)],'-k','linewidth',0.6)
end
p(1) = plot(c2_1(:,2),c2_1(:,3).^2,'ok','markersize',10,'markerfacecolor','b');

expf = @(c,d) c(1).*exp(c(2).*d) + c(3).*exp(c(4).*d) +c(5);

a0 = [87.7984   -0.0200  342.5266   -0.0016  284.2305];
% time constants fixed around 45 and 650 us per Kathuria et all 2014
a0min = [5 a0(2).*1.1 1 a0(4)*1.1 16.9^2]; 
a0max = [inf a0(2).*0.9 inf a0(4).*0.9 17.1^2];
fit1 = lsqcurvefit(expf,a0,c2_1(:,2),c2_1(:,3).^2,a0min,a0max);
xmodel = 0:1:4000;
ymodel1 = expf(fit1,xmodel);
plot(xmodel,ymodel1,'--b','linewidth',2)

c2_2e = c2_2(:,3).^2 .* 2 .* (c2_2(:,4)./c2_2(:,3));
for i = 1:length(c2_2(:,2));
    plot([c2_2(i,2) c2_2(i,2)], [c2_2(i,3).^2+c2_2e(i) c2_2(i,3).^2-c2_2e(i)],'-k','linewidth',0.6)
end
p(2) = plot(c2_2(:,2),c2_2(:,3).^2,'ok','markersize',10,'markerfacecolor','r');

expf = @(c,d) c(1).*exp(c(2).*d) + c(3).*exp(c(4).*d) +c(5);

a0 = [87.7984   -0.0200  342.5266   -0.0016  284.2305];
% time constants fixed around 45 and 650 us per Kathuria et all 2014
a0min = [5 a0(2).*1.1 1 a0(4)*1.1 16.5.^2]; 
a0max = [inf a0(2).*0.9 inf a0(4).*0.9 17.1^2];
fit2 = lsqcurvefit(expf,a0,c2_2(:,2),c2_2(:,3).^2,a0min,a0max);
xmodel = 0:1:4000;
ymodel2 = expf(fit2,xmodel);
plot(xmodel,ymodel2,'--r','linewidth',2)



format_figure

ylabel('R_G^2 (Å^2)')
xlabel('time (seconds)')
box on

legend(p,'mixer 1','mixer 2')
legend box off


%% plot I0

i01 = importdata('08_2019_avg_cytc_avg_I0s.csv')
i01 = i01.data
i02 = importdata('12_2019_avg_cytc2_avg_I0s.csv')
i02 = i02.data

figure
hold on

for i = 1:length(i01(:,2))
   e(i)= plot([i01(i,2) i01(i,2)],[i01(i,3)./i01(1,3)+i01(i,4) i01(i,3)./i01(1,3)-i01(i,4)],'-k')
end

plot(i01(:,2),i01(:,3)./i01(1,3),'ok','markersize',10,'markerfacecolor','r')

for i = 1:length(i02(:,2))
   e(i)= plot([i02(i,2) i02(i,2)],[i02(i,3)./i02(1,3)+i02(i,4) i02(i,3)./i02(1,3)-i02(i,4)],'-k')
end

plot(i02(:,2),i02(:,3)./i02(1,3),'ok','markersize',10,'markerfacecolor','b')

format_figure

box on

xlabel('microseconds')
ylabel('I_0 (normalized)')
