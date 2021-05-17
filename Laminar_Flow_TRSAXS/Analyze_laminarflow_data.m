%script to read the raw data file in each NaCl data directory "data.mat"
%and output:
%   1. normalization
%   2. assembly plot fit to Wiebull CDF
%   3. overlay of data with Wiebull PDF

%% find form factor


clear('agg','aggT','aggA','aggA1','aggA2','dat','datt')
datb = open('data.mat');
datc =datb.datc;

dat = datc{1}
EQ_GC = @(beta,Q)beta(3) + beta(2).*(2.* (exp(-(Q.*beta(1)).^2) + (Q.*beta(1)).^2 - 1)./(Q.*beta(1)).^4);
nm = 2
datt = dat(:,:,nm);
q_max = 0.2;
q_min = 0.02;

[ind m] = find(dat(:,1)>q_min & dat(:,1)<q_max);

q = datt(:,1);
i = datt(:,2);
beta_min = [22 0 0];
beta_max = [35 inf 0.005]
beta0 = [25 max(i) 0.00001];

OPTIONS=optimset;
OPTIONS=optimset(OPTIONS,'HessUpdate','steepdesc');
OPTIONS.MaxIter = 500000;
OPTIONS.maxFunEvals = 100000; 
OPTIONS=optimset(OPTIONS,'maxFunEvals',100000);
OPTIONS=optimset(OPTIONS,'tolX',1e-17);
OPTIONS=optimset(OPTIONS,'tolFun',1e-17);
disp('running')
[betaf, resnorm, residual, exitflag,output,lambda,jacobian]  = lsqcurvefit(EQ_GC,beta0,q(ind),i(ind),beta_min,beta_max,OPTIONS);
%[betaf, resnorm, residual, exitflag,output,lambda,jacobian]  = lsqcurvefit(EQ_GC,beta0,q(ind),i(ind));




q_model= logspace(log10(10E-4),log10(0.5));
i_model= EQ_GC(betaf,q_model);
f1 = figure
hold on
[a b] = logsmooth(datt(:,1),datt(:,2),50);
plot(a,b,'ok','markersize',5,'markerfacecolor',[1 0 0])
p1(2) = plot(q_model,i_model,'--r','linewidth',2);
set(gca,'yscale','log','xscale','log')

%calc assembly

% use the following vc values
% 500 mM 0.055

vc = 0.045;
qmv = ones(1,16)*vc;

q_max = 0.4;

q_t1 = 0.008;
q_t2 = 0.01;
%for 200 mM
q_t1 = 0.01;
q_t2 = 0.012;
dat = 1;
set(0,'currentfigure',f1)
cm = jet(48);
for n1 = 1:16
    clear dat
    dat = datc{n1}
    
    for n2 = 1:length(dat(1,1,:));
        datt = dat(:,:,n2);
        [ind m] = find(datt(:,1)>qmv(n2) & datt(:,1)<q_max);
        q = datt(:,1);
        i = datt(:,2);
        er = datt(:,3);
        betaz = betaf;
        %betaz(2) = 1;
        norm = i(ind)\EQ_GC(betaz,q(ind));
        i0 = betaf(2);
        
        [tind mt] = find(dat(:,1)>q_t1 & dat(:,1)<q_t2);
        v1 = norm * median(i(tind));
        v1e = (sqrt(sum(er(tind).^2))./length(er(tind)))*norm;
        v2 = mean(EQ_GC(betaz,q(tind)));
        aggt = v1 -v2;
        aggz(n2) = aggt;
        agget(n2) = v1e;
        [a b] = logsmooth(datt(:,1),datt(:,2),50);
        if n1 ==5;
        plot(a,b.*norm,'ok','markersize',5,'markerfacecolor',cm(n1*3-1,:))
        end
        
    end
    agg{n1} = aggz;
    agge{n1} = agget;
    aggA(n1) = mean(aggz);
    aggAe(n1) = sqrt(sum(agget.^2))./length(agget);
    
end
set(gca,'linewidth',2,'fontsize',18,'fontweight','bold')
ylabel('I(q)')
xlabel('q (A^{-1})')


% plot agg for fast flow
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
%plot Assembly and CDF

figure
hold on



for z = 1:16
    aggT = agg{z};
    vt = 1:length(aggT);
    %vt = [1 2 3 4 5]
    val = ones(1,length(aggT(vt)));
    
    plot(time(z)*val,aggT(vt),'ok','markersize',4,'markerfacecolor','k')
    plot(time(z),mean(aggT(vt)),'ok','markersize',8,'markerfacecolor','r')
    aggA1(z) = mean(aggT(vt));
    %aggA2(z) = max(aggT) - min(aggT);
    aggA2(z) = std(aggT);
end

box on
format_figure
xlabel('time (ms)')
ylabel('assembly')


%weibull
fun = @(beta,x)beta(1) * (1 - exp(-(x/beta(2)).^beta(3)));
pdf = @(beta,x) (beta(3)/beta(2))*(x./beta(2)).^(beta(3)-1) .* exp(-(x/beta(2)).^beta(3));

beta0 = [max(aggA1) 0.001 30 0];
beta_min = [0 0 0 0];
beta_max = [inf inf inf 0.00000001];

% beta0w = [max(aggA1) b1 b2];
% beta_minw = [0 b1*0.9999 b2*0.9999];
% beta_maxw = [inf b1*1.00001 b2*1.00001];

%parameters from 500mM NaCl fit
bw500 = [0.0087   20.8174    1.7619];
%400mM
bw400 = [ 0.0057   23.1130    2.0317];
%300mM
bw300 = [0.0025   50.5465    1.8369];
%200mM extrapolated parameters
bw200 = [0.002 378.4769 1.8];

beta0w = bw300;
beta_minw = [0 0 1.76];
beta_maxw = [inf 380 inf];


xmodel = 0:0.1:1400;


[betafnw, resnormw, residualw, exitflagw,outputw,lambdaw,jacobianw]  = lsqcurvefit(fun,beta0w,time,aggA1',beta_minw,beta_maxw);

ymodelw = fun(betafnw,xmodel);



plot(xmodel,ymodelw,'--k')

xlim([0 max(time)+1])
vline(betafnw(2))

%plot PDF
fh = figure
lc = [0 0 0];
rc = [1 0 0];
set(fh,'defaultAxesColorOrder',[lc; rc]);
hold on
yyaxis left
plot(xmodel,pdf(betafnw,xmodel),'--r','linewidth',2)
xlim([0 max(time)+100])
ylabel('P')
yyaxis right
plot(time,aggA2,'ok','markersize',8,'markerfacecolor','r')
format_figure
xlabel('time (ms)')
box on
%xlim([0 max(time)+1])
ylim([0 5*max(aggA2)])
%xlim([0 1200])
