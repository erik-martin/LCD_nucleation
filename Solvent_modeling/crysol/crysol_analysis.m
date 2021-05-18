cd files
% This script requires a function to preform Guinier analysis (dis_guinier3)

%low contrast data
h1 = importdata('h102.abs',' ',1);
h1 = h1.data;
h2 = importdata('h202.abs',' ',1);
h2 = h2.data;
m1 = importdata('m102.abs',' ',1);
m1 = m1.data;
m2 = importdata('m202.abs',' ',1);
m2 = m2.data;
l1 = importdata('l102.abs',' ',1);
l1 = l1.data;
l2 = importdata('l202.abs',' ',1);
l2 = l2.data;

%high contrast data
h1b = importdata('h101.abs',' ',1);
h1b = h1b.data;
h2b = importdata('h201.abs',' ',1);
h2b = h2b.data;
m1b = importdata('m101.abs',' ',1);
m1b = m1b.data;
m2b = importdata('m201.abs',' ',1);
m2b = m2b.data;
l1b = importdata('l101.abs',' ',1);
l1b = l1b.data;
l2b = importdata('l201.abs',' ',1);
l2b = l2b.data;

md = importdata('modeldimer00.abs',' ',1);
md = md.data;

q = h1(:,1);
dat =[h1(:,2) h2(:,2) m1(:,2) m2(:,2) l1(:,2) l2(:,2) md(:,2)];
dat2 =[h1b(:,2) h2b(:,2) m1b(:,2) m2b(:,2) l1b(:,2) l2b(:,2)];
dat3 =[h1b(:,2) h2(:,2) m1(:,2) m2(:,2) l1b(:,2) l2b(:,2)];


% figure
% hold on
% cm = jet(6)
% for i = 1:6
%     plot(q,dat(:,i),'-','color',cm(i,:));
% end
% plot(q,dat(:,7),'--k')
% plot(q,mean(dat,2),'-k','linewidth',2)
% 
% set(gca,'yscale','log','xscale','log')
% 
figure
hold on
wt = [1 0.9 0.75 0.5 0.25 0.1 0]

cm = jet(7);
for i = 1:7
    plot(q,wt(i).* mean(dat(:,1:6),2) + (1-wt(i)).* dat(:,7),'-','color',cm(i,:));
    rgio(:,i) = dis_guinier3([q wt(i).* mean(dat(:,1:6),2) + (1-wt(i)).* dat(:,7)],22-i);
    close 
end
set(gca,'yscale','log','xscale','log')
box on
legend('0','10%','25%','50%','75%','100%')
format_figure
ylabel('I(q)')
xlabel('q')

figure
yyaxis left
plot(1-wt,rgio(1,:),'ok','markersize',8,'markerfacecolor','b')
ylabel('R_G')
yyaxis right
plot(1-wt,100*(rgio(2,:)-rgio(2,1))/rgio(2,1),'ok','markersize',8,'markerfacecolor','r')
ylabel('% increase in I_0')

xlabel('fraction dimer')
format_figure
figure
hold on
for i = 1:6
    plot(q*rgio(1,i),(q*rgio(1,i)).^2 .* (wt(i).* mean(dat(:,1:6),2) + (1-wt(i)).* dat(:,7))/rgio(2,i),'-','color',cm(i,:));
end

box on
legend('0','10%','25%','50%','75%','100%')

ylabel('(q R_G)^2 * (I(q)/I_0)')
xlabel('q R_G')

format_figure

figure
hold on

plot(q,mean(dat(:,1:6),2),'-r','linewidth',2)
plot(q,mean(dat2,2),'-b','linewidth',2)
plot(q,mean(dat3(:,5:6),2),'-k','linewidth',2)
rgio2(:,1) = dis_guinier3([q mean(dat(:,1:6),2)],22);
close
rgio2(:,2) = dis_guinier3([q mean(dat2,2)],22);
close
rgio2(:,3) = dis_guinier3([q mean(dat3(:,5:6),2)],22);
close
set(gca,'yscale','log','xscale','log')
box on
format_figure
legend('low contrast','high contrast','high contrast/compact')
ylabel('I(q)')
xlabel('q')

figure
hold on
plot(q*rgio2(1,1),(q*rgio(1,1)).^2 .* mean(dat(:,1:6),2)/rgio(2,1),'-r');
plot(q*rgio2(1,2),(q*rgio(1,2)).^2 .* mean(dat2,2)/rgio(2,2),'-b');
plot(q*rgio2(1,3),(q*rgio(1,3)).^2 .* mean(dat3(:,5:6),2)/rgio(2,3),'-k');
box on
legend('low contrast','high contrast','high contrast/compact')

ylabel('(q R_G)^2 * (I(q)/I_0)')
xlabel('q R_G')

format_figure


figure
yyaxis left
bar([1:3]-0.1,rgio2(1,:),'barwidth',0.2)
ylabel('R_G')
ylim([20 28])
yyaxis right
bar([1:3]+0.1,rgio2(2,:),'barwidth',0.2)
ylabel('I_0')
ylim([0.01 0.02])
set(gca,'xtick',[1 2 3],'xticklabel',{'low contrast','high contrast','high contrast/compact'})
xtickangle(45)
format_figure
