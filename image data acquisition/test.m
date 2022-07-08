close all

P = [p; flipud(p)];
dist_to_center = flipud(w-d);
D = [d; (flipud(d) +2.*dist_to_center)]; % Mirror around lens center
D = D - w; % shift the hwole set to be centered around 0

fp2m = fit(D,P,'poly2','Lower',[-inf -inf,-inf],'Upper',[inf inf,inf]);
xp2m = linspace(-w,w,1000)';
yp2m = feval(fp2m,xp2m);

fp3m = fit(D,P,'poly6','Lower',[-inf -inf,-inf],'Upper',[inf inf,inf]);
xp3m = linspace(-w,w,1000)';
yp3m = feval(fp3m,xp3m);

% fp2 = fit(D(1:sind-1),P(1:sind-1),'poly2','Lower',[-inf -inf,-inf],'Upper',[inf inf,inf]);
% xp2 = linspace(-w,0,1000)';
% yp2 = feval(fp2,xp2);

fp4m = fit(D,P,'poly4','Lower',[-inf -inf,-inf -inf 1.5]);
xp4m = linspace(-w,w,1000)';
yp4m = feval(fp4m,xp4m);


sind = length(p) + 1;
figure(1)
set(gcf,'position',[20,50,1300,1000])
hold on
scatter(D(1:sind-1),P(1:sind-1))
plot(xp2m(1:500),yp2m(1:500),'r')
plot(xp3m(1:500),yp3m(1:500),'k')
plot(xp4m(1:500),yp4m(1:500),'m')
grid minor


