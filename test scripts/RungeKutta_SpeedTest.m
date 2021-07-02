% Testing the speed of Runge kutta solutions on the graded index function
% as well as the Meggit Rochow solver and a simple numerical Snell border
% solver.

% Results MR is ~50% faster than RK
%         SN is ~800% faster than RK

[rk1,rk2,SN,MR] = spedTest(100000);
disp(rk1)
% disp(rk2)
disp(SN)
disp(MR)

disp(strcat('RK1/SN =',num2str(rk1/SN)));
disp(strcat('RK1/MR =',num2str(rk1/MR)));


function [rk1,rk2,SN,MR] = spedTest(nRuns)
n0 = 1.3;
n1 = 1.45;
k = n0-n1;
n = @(r) k*(r(1)^2 + r(2)^2) + n1;
ds = 10^-6;
p0 = [0.7 2];
v0 = [0 -1];
r = 1;
[ip0,~,~] = circleIntersect(r,p0,v0);
ip0 = ip0 + [0 -0.4];
D = @(r) 2*k.*[k.*(r(1)^3 + r(1)*r(2)^2) + n1*r(1);
               k.*(r(2)^3 + r(2)*r(1)^2) + n1*r(2)];
           dnds = @(r) 2*k*[p0(1);p0(2)];
           D2 = @(r) n(r).*dnds(r);
tic
for i = 1:nRuns
    A = ds.*D2(p0);
    B = ds.*D2(p0 + ds/2.*v0 + ds/8.*A);
    C = ds.*D2(p0 + ds.*v0 + ds/2.*B);
    p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
    v1 = v0 + 1/6.*(A + 4.*B + C);
end
rk1 = toc;
rk2 = 1;
% tic
% for i = 1:nRuns
%     A = ds.*D(p0);
%     B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
%     C = ds.*D(p0 + ds.*v0 + ds/2.*B);
%     p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
%     v1 = v0 + 1/6.*(A + 4.*B + C);
% end
% rk2 = toc;

k = n0-n1;
nFunc = @(r) k*r.^2 + n1;

tic
for i = 1:nRuns
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
SN = toc;

% v = v./norm(v);
k = n0-n1;
nFunc = @(r) k*r^2 + n1;
r0 = 1;
tic
for i = 1:nRuns
    sN = -p0';
    sN = sN./sqrt(sN(1)^2 + sN(2)^2);
    n = nFunc(r0);
    theta = acos(sN(1)*v0(1) + sN(2)*v0(2));
    dn = -2*k*r0;
    r = (n) / (sin(theta)*dn);
    [p1,v1] = curvePath(p0,v0,r,ds);
end
MR = toc;


function [p1,v1] = rungekuttaTrace1(p0,v0)
% dnds = @(r) 2*k*[p0(1);p0(2)];
% D2 = @(r) n(r).*dnds(r);
A = ds.*D2(p0);
B = ds.*D2(p0 + ds/2.*v0 + ds/8.*A);
C = ds.*D2(p0 + ds.*v0 + ds/2.*B);            
p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
v1 = v0 + 1/6.*(A + 4.*B + C);
end

function [p1,v1] = rungekuttaTrace2(p0,v0)       
A = ds.*D(p0);
B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
C = ds.*D(p0 + ds.*v0 + ds/2.*B);
p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
v1 = v0 + 1/6.*(A + 4.*B + C);
end

end


