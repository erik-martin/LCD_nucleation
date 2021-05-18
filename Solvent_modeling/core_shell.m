%% run
% core shell


%beta1 = core radius
%beta2 = shell thickness
%beta3 = core sld
%beta4 = shell sld
%beta5 = solvent sld

% 4 * 3.14 * (beta(1) + beta(2))^3 * 

FQ = @(beta,q) 9/(4 * 3.14 * (beta(1) + beta(2))^3) * ( ... 
    ( ( 3.14 * 4 / 3  * beta(1)^3 ) ...
    * (beta(3) - beta(4)) ...
    * (sin(q .* beta(1)) - q .* beta(1) .* cos(q * beta(1))) ./ (q * beta(1)).^3 ) ...
    + ( ( 3.14 * 4 / 3  * (beta(1) + beta(2))^3 ) ...
    * (beta(4) - beta(5)) ...
    * (sin(q * (beta(1) + beta(2))) - q .* (beta(1) + beta(2)) .* cos(q * (beta(1) + beta(2)))) ./ (q * (beta(1) + beta(2))).^3 ) ...
    );


q = 0.005:0.001:0.3;
%beta1 = [20 15 4 1.1 1];
beta1 = [20 10 1E-6 2E-6 1E-6];
figure
hold on
cm = jet(21);
for i = 1:5;
    %modify core density
    
    beta1(3) = beta1(4)*(1 + (i-1)*0.1);
    beta1(3)/beta1(4)
    plot(q,3/(4* 3.14 * (beta1(1) + beta1(2))^3) * FQ(beta1,q).^2,'-','color',cm(5*i-4,:))
    
    %the following command requires calculating the Rg/I0 using Guinier
    %analysis
    %rgio(:,i)= dis_guinier3([q' (3/(4* 3.14 * (beta1(1) + beta1(2))^3) * FQ(beta1,q).^2)'],20)
    %close
end

set(gca,'yscale','log','xscale','log')
xlim([0.008 0.35])
%legend toggle
box on
format_figure
ylabel('I(q)')
xlabel('q')
colormap(cm)
%colorbar


colorbar('location','northoutside','Ticks',[0.1,0.9],...
         'TickLabels',{'small','large'})
     
% If the Guinier Rg/I0 is calculated, this can be uncommented to make panel B 

% figure
% yyaxis left
% plot([1:5]*0.1,rgio(1,:),'ok','markersize',8,'markerfacecolor','b')
% ylabel('R_G')
% yyaxis right
% plot([1:5] *0.1,100*(rgio(2,:)-rgio(2,1))/rgio(2,1),'ok','markersize',8,'markerfacecolor','r')
% ylabel('% increase in I_0')
% 
% xlabel('excess core density')
% format_figure
% 