
%% run
% make 3D plot of raw data 
call

datc = open('data.mat');
datc =datc.datc;

sz = size(datc);

for i = 1:sz(2);
    tmp = datc{i};
    %randomly chose the second file
    matc(:,:,i) = tmp(:,:,2);
    Z(:,i) = tmp(:,2,1); 
end

figure
hold on
% [a b] = meshgrid(matc(:,1,1),1:16);
% mesh(a,b,Z','o')
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
cm= jet(5*length(time));
for i = 1:sz(2);
    [x1 y1] = logsmooth(matc(:,1,i),matc(:,2,i),40);
    plot3(matc(:,1,i),time(i)*ones(size(matc(:,2,i))),matc(:,2,i),'-','linewidth',0.5,'color',cm(i*5-4,:))
    plot3(x1,time(i)*ones(size(y1)),y1,'ok','markerfacecolor',cm(i*5-4,:))
    edv(:,i) = [y1(end) y1(1)]; 
end
set(gca,'zscale','log','xscale','log')
view([10 10])
zlim([median(edv(1,:)) max(edv(2,:))])
xlim([0 0.3])
%box on
format_figure
ylabel('time (ms)')
xlabel('q (Å^{-1})')
zlabel('I(q)')

%% run
% make overlay of single data point
call

datc = open('data.mat');
datc =datc.datc;

pt = 2;
datd = datc{pt};
sz = size(datd);
cm = jet(sz(3)*4);

figure
hold on
for i = 1:sz(3);
    [x1 y1] = logsmooth(datd(:,1,i),datd(:,2,i),40);
    plot(datd(:,1,i),datd(:,2,i),'-','linewidth',0.5,'color',cm(i*4-3,:))
    plot(x1,y1,'ok','markerfacecolor',cm(i*4-3,:))
    edv(:,i) = [y1(end) y1(1)];
end
    
set(gca,'yscale','log','xscale','log')
xlim([0 0.3])
ylim([median(edv(1,:)) max(edv(2,:))])
box on
format_figure

xlabel('q (Å^{-1})')
ylabel('I(q)')