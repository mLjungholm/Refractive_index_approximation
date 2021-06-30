% x = linspace(0,1,1000);
% tline1 = grin.dX(1000,1000:end);
% xline = grin.X(1000,1000:end);
% tline2 = 2.*k.*x;

% close all
% 
% figure(1)
% hold on; grid on
% plot(x,tline2+0.001,'r')
% plot(xline,tline1,'b')
% 
% nline = grin.P(1000,1000:end);
% figure(2)
% hold on; grid on
% plot(x,nGradient(x)+0.001,'r')
% plot(xline,nline,'b')

n0 = 1.3;
n1 = 1.45;
k = n0-n1;
r = 1;
ds = 10^-5;
nFunc = @(x) k.*x.^2 + n1;
p0 = [0.7 2];
v0 = [0 -1];

[ip0,~,~] = circleIntersect(r,p0,v0);
ip0 = ip0 + [0 -0.4];
r0 = norm(ip0);

tic
dr = getLocalRadiusKnownGradient(ip0,v0,r0,n0,n1);
[p1,v1] = curvePath(ip0,v0,dr,ds);
toc
tic
[p2,v2] = rungekuttaTrace(ip0',v0',ds,n0,n1);
toc

dang = @(u,v) acos(u(1)*v(1) + u(2)*v(2)/(sqrt(u(1)^2 + u(2)^2)*sqrt(v(1)^2 + v(2)^2)));

ang1 = dang(v0,v1);
ang2 = dang(v0,v2);
difAng = ((ang1-ang2)/ang1).*100;

function r = getLocalRadiusKnownGradient(p,v,r0,n0,n1)
sN = -p';
sN = sN./sqrt(sN(1)^2 + sN(2)^2);
v = v./norm(v);
k = n0-n1;
nFunc = @(r) k*r^2 + n1;
n = nFunc(r0);
theta = acos(sN(1)*v(1) + sN(2)*v(2));
dn = -2*k*r0;
r = (n) / (sin(theta)*dn);
end

function [p1,v1] = rungekuttaTrace(p0,v0,ds,n0,n1)
v0 = v0./norm(v0);
k = n0-n1;
n = @(r) k*(r(1)^2 + r(2)^2) + n1;
dnds = @(r) 2*k*[p0(1);p0(2)];
D2 = @(r) n(r).*dnds(r);
D = @(r) 2*k.*[k.*(r(1)^3 + r(1)*r(2)^2) + n1*r(1);
               k.*(r(2)^3 + r(2)*r(1)^2) + n1*r(2)];

A = ds.*D2(p0);
B = ds.*D2(p0 + ds/2.*v0 + ds/8.*A);
C = ds.*D2(p0 + ds.*v0 + ds/2.*B);            
% A = ds.*D(p0);
% B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
% C = ds.*D(p0 + ds.*v0 + ds/2.*B);
p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
v1 = v0 + 1/6.*(A + 4.*B + C);
end

