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
ds = 10^-6;
nFunc = @(x) k.*x.^2 + n1;
p0 = [0.7 2];
v0 = [0 -1];

[ip0,~,~] = circleIntersect(r,p0,v0);
[dr,nn1,dn1] = getLocalRadiusKnownGradient(ip0,v0,r,nFunc,k);
[p1,v1] = curvePath(ip0,v0,dr,ds);
[p2,v2,nn2,dn2] = rungekuttaTrace(ip0',v0',ds,n0,n1);

dang = @(v0,v1) acos(v0(1)*v1(1) + v0(2)*v1(2)/(sqrt(v0(1)^2 + v0(2)^2)*sqrt(v1(1)^2 + v1(2)^2)));

ang1 = dang(v0,v1);
ang2 = dang(v0,v2);

difAng = ((ang1-ang2)/ang1).*100;

function [r,n,dn] = getLocalRadiusKnownGradient(p,v,r0,nFunc,k)
sN = -p';
sN = sN./sqrt(sN(1)^2 + sN(2)^2);
n = nFunc(r0);
theta = acos(sN(1)*v(1) + sN(2)*v(2));
dn = -2*k*r0;
r = (n) / (sin(theta)*dn);
end

function [p1,v1,np,dn] = rungekuttaTrace(p0,v0,ds,n0,n1)
v0 = v0./norm(v0);
k = n0-n1;
n = @(r) k*sqrt(r(1)^2 + r(2)^2)^2 + n1;
dnds = @(r) r.*[k/sqrt(r(1)^2 + r(2)^2) ; k/sqrt(r(1)^2 + r(2)^2)];
D = @(r) n(r).*dnds(r);
A = ds.*D(p0);
B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
C = ds.*D(p0 + ds.*v0 + ds/2.*B);
p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
v1 = v0 + 1/6.*(A + 4.*B + C);

np = n(p0);
dn = norm(dnds(p0));
end

