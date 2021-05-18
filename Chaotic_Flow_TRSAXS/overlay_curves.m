%script to overlay and normalize SAXS curves

%run from the folder 'dats'
cd A1LCD/dats
clear
files = dir('*.dat');
num = length(files);
for i = 1:num
    dtmp = importdata(files(i).name,' ',4);
    dat(:,:,i) = dtmp.data;
    int(i) = trapz(dat(:,1,i),dat(:,2,i));
    [a b] = logsmooth(dat(:,1,i),dat(:,2,i),30);
    dats(:,:,i) = [a' b'];
end

i0 = dlmread('../mff/I0.dat')
rg = dlmread('../mff/Rg.dat')
% plot
skip = 2;
cm = jet(num);
figure
hold on
for i = 1:skip:num
    plot(dats(:,1,i).*rg(i),dats(:,2,i)./i0(i),'-o','color',cm(i,:),'markerfacecolor',cm(i,:))
end
format_figure
ylabel('I(q)/I_0')
xlabel('q * R_G')
set(gca,'yscale','log','xscale','log')
box on

