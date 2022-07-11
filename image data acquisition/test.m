close all

P = [p; flipud(p)];
dist_to_center = flipud(w-d);
D = [d; (flipud(d) +2.*dist_to_center)]; % Mirror around lens center
D = D - w; % shift the hwole set to be centered around 0

% fp2m = fit(D,P,'poly2','Lower',[-inf -inf,-inf],'Upper',[inf inf,inf]);
% xp2m = linspace(-w,w,1000)';
% yp2m = feval(fp2m,xp2m);

fp6m = fit(D,P,'poly6');
xp6m = linspace(-w,w,1000)';
yp6m = feval(fp6m,xp6m);

sind = length(p) + 1;
phaseD = -D(1:sind-1);
phaseS = P((1:sind-1));

newEdge = xp6m(find(yp6m > 0,1,'last'));

figure(1)
set(gcf,'position',[20,50,1300,1000])
hold on
scatter(phaseD,phaseS)
plot(xp6m(500:end),yp6m(500:end),'k')
plot([0;0],[0;max(yp6m)],'k')
plot([0;max(xp6m)],[0;0],'k')
plot(newEdge,0,'o')
grid minor


