

%default parameters
pars(1) = 600; %radius of aggregate
pars(2) = pars(1) *0.1; % fuzziness 
pars(3) = 0.00005; % vf of spheres
pars(4) = 0.15; % sphere contrast
pars(5) = 28; % radius of monomer
pars(6) = 0.001; % vf of polymer
pars(7) = 0.15; % polymer contrast
pars(8) = .3 % polydispersity
pars(9) = pars(1) % center of polydisperse R




% fuzzy sphere
% https://pubs.acs.org/doi/pdf/10.1021/la049518x
% beta1 = Radius of sphere measured to where density drops to 1/2
% beta2 = width of the fuzzy region

A1 = @(beta,q) (( 3.*( sin( beta(1).*q ) - beta(1).*q.*cos( beta(1).* q ) ))./ (beta(1).*q).^3) .* exp(-(beta(2) * q)/2);

%polydispersity
% beta8 = polydispersity


p1 = @(beta,r) 1/(sqrt(2*3.14*beta(8).^2*beta(1)^2)) * exp(-(r - beta(9)).^2/(2*beta(8).^2*beta(1).^2));


% fuzzy sphere I0
% beta3 = volume fraction of fuzzy spheres
% beta4 = contrast

I1 = @(beta,q) beta(3) * (3.14 * 4 / 3 * beta(1)^3) * beta(4)^2 ;

P1 = @(beta,q) I1(beta,q) * A1(beta,q).^2;

% gaussian coil
% beta5 = polymer Rg
%

A2 = @(beta,q)(2.* (exp(-(q.*beta(5)).^2) + (q.*beta(5)).^2 - 1)./(q.*beta(5)).^4);

% Gaussian coil I0
% beta6 = polymer volume fraction
% beta7 = polymer contrast


%molecular weight
Ma = 12500;
Na = 6.02E23;
dnp = 1.35E-24; %g/A^3
Vo = Ma / (Na * dnp);
Vob = Ma *2 / (Na * dnp);

I2 = @(beta,q) beta(6) * Vo * beta(7)^2;
I2b = @(beta,q) beta(6) * Vob * beta(7)^2;

P2 = @(beta,q) I2(beta,q) * A2(beta,q);
P2b = @(beta,q) I2b(beta,q) * A2(beta,q);


q = logspace(log10(0.008),log10(0.3),200);
cv1 = P1(pars,q);
cv2 = P2(pars,q);
cvf = P1(pars,q) + P2(pars,q);
% figure
% hold on
% plot(q,cvf,'ok','markerfacecolor','k')
% %plot(q,cv1,'-b','markerfacecolor','b')
% plot(q,cv2,'-r','markerfacecolor','r')
%set(gca,'yscale','log','xscale','log')


r = 0.1:10*pars(1);
cv1p = zeros(1,length(q));
for i = 1:length(r)
    pars(1) = r(i);
    cv1p = cv1p + p1(pars,r(i)) * P1(pars,q);
end


ratio = [0 0.01 0.1 0.25 0.5 1];
cm = jet(length(ratio)*5);


q_t1 = 0.008;
q_t2 = 0.01;
[m ind] = find(q>q_t1 & q<q_t2);
figure
hold on

for i = 1:length(ratio)
    % assume that the maximum vf ratio is 0.02 based on binodal
    val2 = 0.02*ratio(i);
    vf1 = (1-val2);
    vf2 = pars(6) * val2;
    
    %include 1/pars(3) to correct for FC already having a vf
    
    cvfp = P2(pars,q) *vf1 + cv1p*vf2/pars(3);
    plot(q,cvfp,'-','color',cm(i*5-4,:))
    agg(i) = mean(cvfp(ind)) - mean(P2(pars,q(ind)));
    
    
end

set(gca,'yscale','log','xscale','log')
ylim([min(cv2) 2*max(cvfp)])
format_figure
box on
xlabel('q (Å^{-1})')
ylabel('I(q)')

figure
hold on
[a b] = simple_linear_fit(0.02*ratio,agg);
plot(b(:,1),b(:,2),'--k')
plot(0.02*ratio,agg,'ok','markersize',8,'markerfacecolor','r')

format_figure
box on
xlabel('\phi_{aggregate} / \phi_{monomer}')
ylabel('assembly metric')




