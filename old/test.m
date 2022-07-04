close all
p = l.PP(l.Pinc);
d = l.PD(l.Pinc);
w = l.centerLine - l.leftEdge;

% Mirror the data at the lens center
P = [p; flipud(p)];
dist_to_center = flipud(w-d);
D = [d; (flipud(d) +2.*dist_to_center)];
D = D - w;

restr = [1.7 1.75];
f = fit(D,P,'poly2','Lower',[-inf -inf,restr(1)],'Upper',[inf inf,restr(2)]);
% f = fit(D,P,'poly2');
x = linspace(-w,w,1000)';

% x = linspace(0,2*w,1000)';
y = feval(f,x);
miny = min(y(y>0));

% xpos = find(y == miny);
xpos = find(y(x>0) == miny) + find(x>0,1,'first') - 1;
ycheck = miny == y(xpos);
newedge = x(xpos);

x2 = linspace(0,newedge,1000);
y2 = feval(f,x2);

figure(1)
hold on
% scatter(d,p)
% scatter(D,P)
% plot(x,y,'r')
plot(x2,y2)
scatter(D(length(p)+1:end),P(length(p)+1:end))
% plot(x(500:end),y(500:end),'r')
% scatter(D(1:length(p)),P(1:length(p)))
% plot(x(1:500),y(1:500),'r')
grid minor






% uval = -inf;
% f2 = fit(d,p,'power2','Lower',[-inf -inf,uval],'Upper',[0 inf inf]);
% y2 = feval(f2,x(2:end));

% plot(x(2:500),y2(2:500),'g')
% scatter(D,P)
% plot(x,y)
