function [ rgio,qs,ys,fit_dat ] = dis_guinier(varargin)
%UNTITLED7 Summary of this function goes here
%   Function that takes an estimate of the Rg as an input and returns
%   Guinier analysis 
%% demo
if nargin >3
    in_mat = varargin{1};
    reg = varargin{2};
    skip = varargin{3};
    qrg = varargin{4};
    
elseif nargin == 3;
    in_mat = varargin{1};
    reg = varargin{2};
    skip = varargin{3};
    qrg = 1;
elseif nargin == 2;
    in_mat = varargin{1};
    reg = varargin{2};
    skip = 0;
    qrg = 1;
elseif nargin == 1;
    in_mat = varargin{1};
    reg = 30;
    skip = 0;
    qrg = 1;
else
    disp('come on man, you gotta give me something...')
end
   
strt = skip + 1;
Q = in_mat(strt:end,1);
Q2 = in_mat(:,1);
I = in_mat(strt:end,2);
I2 = in_mat(:,2);

%find the appropriate q region using the estimate of Rg

q_index = find(Q<(qrg/reg)); 
q_index2 = find(Q2<(qrg*1.2/reg)); 

rem_cycle = [];
z_ref = 1;
cycle = 0;
val = 1;
%while loop removes outliers in the low q. It tests to see if the residual is outside of outlier*std(residuals)  
outlier =1.5;
while z_ref > 0
    cycle = cycle + 1;
    z_ref = 0;
    q_gun = Q(q_index(val:end));
    i_gun = I(q_index(val:end));

    b0 = [ones(length(q_gun),1) q_gun.^2];
    bf = q_gun.^2;
    y0 = log(i_gun);
    %g0 = b0\y0;
    fun = fittype( @(cof,b,x) cof*x +b);
    fit_temp = fit(bf,y0,fun,'StartPoint',[1 1],'Robust','on');
    g0 = [fit_temp.b;fit_temp.cof];
    y1 = b0*g0;
    R = y0 -y1;
%    z_ref = abs(R(1)) > (mean(abs(R)) + outlier*std(R))
%    rem_cycle = [rem_cycle 1];
%    abs(R);
%    (mean(abs(R)) + outlier*std(R))
    ref = abs(R) > (mean(abs(R)) + outlier*std(R));
%     if sum(ref) > 0 && cycle ==1;
%         val = val+max(find(ref==1))+1;
%     
%         z_ref = 1;
%     end
end
%rem_cycle
SEs = sqrt(sum(R.^2)/(length(q_gun)-2))/sqrt(sum((q_gun.^2 - mean(q_gun.^2)).^2));
%t-value
c_int = 95;
int = 1-(1-c_int/100)/2;
crit = tinv(int,(length(q_gun)-2));

ME = crit*SEs;
g1 = g0 + [0;ME];
off1 = [1 median(q_gun.^2)]*g0 - [1 median(q_gun.^2)]*g1;
g2 = g0 - [0;ME];
off2 = [1 median(q_gun.^2)]*g0 - [1 median(q_gun.^2)]*g2;

Rg0 = sqrt(-[g0(2) g1(2) g2(2)]*3);
Rg_err = abs((Rg0(2)-Rg0(3))/2);
disp(sprintf('R_g = %.1f, I_0 = %.2f',Rg0(1),exp(g0(1))))
disp(sprintf('R_g * q = %.3f',Rg0(1)*max(q_gun)))
disp(sprintf('used q range  %f to %f',min(q_gun),max(q_gun)))
if Rg0(1)*max(q_gun) > 1.1
    disp('WARNING, q * R_g too large, use larger R_g estimate')
elseif Rg0(1)*max(q_gun) < 0.9 
    disp('WARNING, q * R_g too small, change R_g estimate')
else
    disp('all good')
end
rgio = [Rg0(1);exp(g0(1))];
fit_dat = y1;
ys = y0;
qs = q_gun;
fh =figure;
s1 = subplot(3,1,[1 2]);
hold on

pl(1) = plot(Q2(q_index2).^2,log(I2(q_index2)),'o','color',[0.4 0.4 0.4],'markersize',6,'linewidth',0.5);
pl(2) = plot(q_gun.^2,y0,'ok','markersize',6,'markerfacecolor','k');

pl(3) = plot(q_gun.^2,b0*g0,'-k','linewidth',2);
pl(4) = plot(q_gun.^2,off1 + b0*g1,'-b','linewidth',2);
pl(5) = plot(q_gun.^2,off2 + b0*g2,'-r','linewidth',2);
box on
xrange = [min(q_gun.^2) max(q_gun.^2)] + [-0.1*range(q_gun.^2), 0.1*range(q_gun.^2)];
yrange0 = [min(y0) max(y0)] + [-0.5*range(y0), 0.5*range(y0)];
xlim(xrange)
ylim(yrange0)
set(gca,'linewidth',2,'fontsize',18,'fontweight','bold','xticklabel',' ')
ylabel('Log(I)')
tx = xrange(1) + 0.6*(range(xrange));
ty = yrange0(1) + 0.75*(range(yrange0));
t1 = text(tx,ty,sprintf('R_G = %0.2f (+/- %0.2f)',Rg0(1),Rg_err),'fontsize',14,'fontweight','bold');
%residuals
s2 = subplot(3,1,3);
hold on

tst = size(in_mat);
if tst(2)>2;
    R = R./in_mat(q_index(val:end),3);
else
    R = R./sqrt(abs(in_mat(q_index(val:end),2)));
end
pl(5) = plot(q_gun.^2,R,'ok','markersize',6,'markerfacecolor','k','linewidth',0.5);
pl(6) = plot(xrange,[0 0],'--k','linewidth',2);
box on
set(gca,'linewidth',2,'fontsize',18,'fontweight','bold')
yrange1 = [min(R) max(R)] + [-0.15*range(R), 0.15*range(R)];
ylim(yrange1)
xlim(xrange)
ylabel('\Delta I / I^{1/2}')
if tst(2)>2;
    ylabel('\Delta I/\sigma')
end

xlabel('q^2')


%kratky(in_mat,rgio)
end

