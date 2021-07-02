p0 = [0.4,0.4];
v0 = [0,-1];
v0 = v0./norm(v0);
% [p0,~] = circleLineIntersect(1,p0',v0');   
[p1,v1] = snellA(p0,v0);
dang = @(u,v) acos(u(1)*v(1) + u(2)*v(2)/(sqrt(u(1)^2 + u(2)^2)*sqrt(v(1)^2 + v(2)^2)));


figure(1)
hold on; axis equal; grid on
plotCircle(1,2*pi,1)
quiver(p0(1),p0(2),v0(1),v0(2))
quiver(p1(1),p1(2),v1(1),v1(2))

function [p1,v1] = snellA(p0,v0)
ds = 10^-7;
n0 = 1.3;
n1 = 1.45;
k = n0-n1;
nFunc = @(r) k*r.^2 + n1;
N = p0;
r0 = norm(p0);
na = nFunc(r0 - ds/2*r0);
nb = nFunc(r0 + ds/2*r0);
a = v0(1)*N(1) + v0(2)*N(2);
if a > 0
    N = -N; 
end
r = na/nb;
v1 = r*v0 + (r*a - sqrt(1- r^2*(1-a^2)))*N;
v1 = v1./norm(v1);
p1 = p0 + ds.*v1;
end

function [p1,v1] = snellB(p0,v0)
N = p0;
r0 = norm(p0);
na = nFunc(r0 - ds/2*r0);
nb = nFunc(r0 + ds/2*r0);
a = v0(1)*N(1) + v0(2)*N(2);
% a = dot(v0,N);
if a > 0
    N = -N;
    
end

% c = -1*dot(N,v0);
r = na/nb;

v1 = r*v0 + (r*a - sqrt(1- r^2*(1-a^2)))*N;

% if isreal(vNew)  % If vNew is imagenary then total internal reflection occured. 
% %     vRefracted = vNew;
%     reflected = false;
% else
%     vNew = v - 2*dot(v,N).*N;   % Code for reflection.
%     vRefracted = vNew;
%     reflected = true;
% end
% [v1, ~] = Snell(v0, N, na, nb);
p1 = p0 + ds.*v1;
end

function [p,intersect] = circleLineIntersect(r,p0,v)
L = -p0;
tca = v(1)*L(1) + v(2)*L(2);
if tca < 0
    intersect = 0;
    p = nan;
    return
end
d = sqrt(L(1)^2 + L(2)^2 - tca^2);
if d < 0 || d > r
    intersect = 0;
    p = nan;
    return
end
thc = sqrt(r^2-d^2);
t0 = tca - thc;
p = p0 + t0*v;
intersect = 1;
end
