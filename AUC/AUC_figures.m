%% run
%plot AUC data

A1 = readmatrix('AUCcurves2.xlsx','Sheet','cofs_measured','Range','A2:B101');
A2 = readmatrix('AUCcurves2.xlsx','Sheet','cofs_measured','Range','D2:E51');
A3 = readmatrix('AUCcurves2.xlsx','Sheet','cofs_measured','Range','G2:H101');
A4 = readmatrix('AUCcurves2.xlsx','Sheet','cofs_measured','Range','J2:K101');
A5 = readmatrix('AUCcurves2.xlsx','Sheet','cofs_measured','Range','M2:N101');
A6 = readmatrix('AUCcurves2.xlsx','Sheet','cofs_measured','Range','P2:Q51');
sz=[10 10 700 250];
figure('pos',sz)
hold on
cm = winter(4);

a(1) = area(A2(5:end,1),A2(5:end,2)./trapz(A2(5:end,1),A2(5:end,2)),'facecolor',cm(1,:));
p(1) = plot(A2(5:end,1),A2(5:end,2)./trapz(A2(5:end,1),A2(5:end,2)),'-','linewidth',2,'color',cm(1,:))

a(2) = area(A5(10:end,1),A5(10:end,2)./trapz(A5(10:end,1),A5(10:end,2)),'facecolor',cm(2,:));
p(2) = plot(A5(10:end,1),A5(10:end,2)./trapz(A5(10:end,1),A5(10:end,2)),'-','linewidth',2,'color',cm(2,:))

a(3) = area(A3(10:end,1),A3(10:end,2)./trapz(A3(10:end,1),A3(10:end,2)),'facecolor',cm(3,:));
p(3) = plot(A3(10:end,1),A3(10:end,2)./trapz(A3(10:end,1),A3(10:end,2)),'-','linewidth',2,'color',cm(3,:))

a(4) = area(A4(10:end,1),A4(10:end,2)./trapz(A4(10:end,1),A4(10:end,2)),'facecolor',cm(4,:));
p(4) = plot(A4(10:end,1),A4(10:end,2)./trapz(A4(10:end,1),A4(10:end,2)),'-','linewidth',2,'color',cm(4,:))

p(5) = plot(A1(10:end,1),A1(10:end,2)./trapz(A1(10:end,1),A1(10:end,2)),':','color',cm(1,:),'linewidth',2);
p(6) = plot(A6(10:end,1),A6(10:end,2)./trapz(A6(5:end,1),A6(5:end,2)),':','color',cm(3,:),'linewidth',2);

set(a,'facealpha',0.2);

set(gca,'fontsize',18,'fontweight','bold')
box on
xlim([0.7 3.5])
legend(p,...
    '50 mM NaCl, 0.3 mg/mL A1-LCD',...
    '50 mM NaCl, 2.3 mg/mL A1-LCD',...
    '50 mM NaCl, 5.6 mg/mL A1-LCD',...
    '50 mM NaCl, 10.1 mg/mL A1-LCD',...
    '0 mM NaCl, 0.35 mg/mL A1-LCD',...
    '0 mM NaCl, 5.6 mg/mL A1-LCD') 
legend box off
ylabel('c_{NI}(s_0) (fringes/S)')
xlabel('sedimentation coefficient S')

%% run 
% plot the Glibert Theory simulations

G1 = readmatrix('AUCcurves2.xlsx','Sheet','isodesmicGilbert','Range','B1:C54');
G2 = readmatrix('AUCcurves2.xlsx','Sheet','isodesmicGilbert','Range','F1:G54');
G3 = readmatrix('AUCcurves2.xlsx','Sheet','isodesmicGilbert','Range','J1:K54');
G4 = readmatrix('AUCcurves2.xlsx','Sheet','isodesmicGilbert','Range','N1:O54');


sz=[10 10 700 250];
figure('pos',sz)
hold on

a1(4) = area(G4(:,1),G4(:,2)./max(G4(:,2)),'facecolor',cm(4,:))
plot(G4(:,1),G4(:,2)./max(G4(:,2)),'color',cm(4,:))

a1(3) = area(G3(:,1),G3(:,2)./max(G3(:,2)),'facecolor',cm(3,:))
plot(G3(:,1),G3(:,2)./max(G3(:,2)),'color',cm(3,:))

a1(2) = area(G2(:,1),G2(:,2)./max(G2(:,2)),'facecolor',cm(2,:))
plot(G2(:,1),G2(:,2)./max(G2(:,2)),'color',cm(2,:))

a1(1) = area(G1(:,1),G1(:,2)./max(G1(:,2)),'facecolor',cm(1,:))
plot(G1(:,1),G1(:,2)./max(G1(:,2)),'color',cm(1,:))


set(a1,'facealpha',0.4);
ylabel('c_{NI}(s_0) (fringes/S)')
xlabel('sedimentation coefficient S')

set(gca,'fontsize',18,'fontweight','bold')
box on
xlim([1 3])
ylim([0 1.1])



G21 = readmatrix('AUCcurves2.xlsx','Sheet','2step-isodesmicGilbert','Range','B1:C54');
G22 = readmatrix('AUCcurves2.xlsx','Sheet','2step-isodesmicGilbert','Range','F1:G54');
G23 = readmatrix('AUCcurves2.xlsx','Sheet','2step-isodesmicGilbert','Range','J1:K54');
G24 = readmatrix('AUCcurves2.xlsx','Sheet','2step-isodesmicGilbert','Range','N1:O54');


figure('pos',sz);
hold on
a2(4) = area(G24(:,1),G24(:,2)./max(G24(:,2)),'facecolor',cm(4,:))
plot(G24(:,1),G24(:,2)./max(G24(:,2)),'color',cm(4,:))

a2(3) = area(G23(:,1),G23(:,2)./max(G23(:,2)),'facecolor',cm(3,:))
plot(G23(:,1),G23(:,2)./max(G23(:,2)),'color',cm(3,:))

a2(2) = area(G22(:,1),G22(:,2)./max(G22(:,2)),'facecolor',cm(2,:))
plot(G22(:,1),G22(:,2)./max(G22(:,2)),'color',cm(2,:))

a2(1) = area(G21(:,1),G21(:,2)./max(G21(:,2)),'facecolor',cm(1,:))
plot(G21(:,1),G21(:,2)./max(G21(:,2)),'color',cm(1,:))





set(a2,'facealpha',0.4);
ylabel('c_{NI}(s_0) (fringes/S)')
xlabel('sedimentation coefficient S')

set(gca,'fontsize',18,'fontweight','bold')
box on
xlim([1 3])
ylim([0 1.1])


%% fits
G31 = readmatrix('AUCcurves2.xlsx','Sheet','sw_isotherm_fits','Range','A3:B6');
G32 = readmatrix('AUCcurves2.xlsx','Sheet','sw_isotherm_fits','Range','D3:E52');
G33 = readmatrix('AUCcurves2.xlsx','Sheet','sw_isotherm_fits','Range','G3:H52');

sz=[10 10 250 250];
figure('pos',sz);
hold on
cm = winter(4);

for i = 1:length(G31(:,1))
    c1(i) = plot(G31(i,1),G31(i,2),'ok','markerfacecolor',cm(i,:),'markersize',10)
end

%d1 = plot(G32(:,1),G32(:,2),'-','color',[0.2 0.2 0.2])
d1 = plot(G33(:,1),G33(:,2),'-','color',[0.2 0.2 0.2])

yl = [0.9*min(G31(:,2)) 1.1*max(G31(:,2))];
xl = [0.5*min(G31(:,1)) 1.5*max(G31(:,1))];
xlim(xl)
ylim(yl)

set(gca,'xscale','log',...
'fontweight','bold',...
'fontsize',18,...
'fontname','helvetica')
box on
ylabel('sw (S)')
xlabel('[A1-LCD] (?M)')


