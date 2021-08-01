% Functions for testing the different solvers for different index profiles

function [p1,v1,alpha,dserr,np,d] = meggit_step_elliptical(ds,p0,v0,n0,n1,ra,rb,nProfile)
switch nProfile
    case 'parabolic'
%         k = (n0-n1)/rn^2;
        n = @(r,rn) (n0-n1)/rn^2*r^2 + n1;
        dndr = @(r,rn) 2*(n0-n1)/rn^2*r;
    case 'linear'
%         k = (n0-n1)/rn;
        n = @(r,rn) (n0-n1)/rn*r + n1;
        dndr = @(r,rn) (n0-n1)/rn;
    case 'eliptical'
        k = (n1-n0)/rn;
        n = @(r) k*sqrt(rn^2 - r^2) + n0;
        dndr = @(r) -r*k/sqrt(rn^2 - r^2);
    case 'step'
        k = n1-n0;
        n = @(r) n0 + k*(1-sign(r-rn)/2);
        dndr = @(r) 0;
    case 'luneburg'
        R = rn;
        n = @(r) sqrt(2-r^2/R^2);
        dndr = @(r) -r/sqrt(2-r^2/R^2)/R^2;
end
if norm(p0) == 0
    testflag = 1;
end

r0 = norm(p0);
gV = -p0./r0;
gV = gV./norm(gV);
theta = acos(dot(v0,gV)./norm(v0));
psi = acos(dot(p0,[1 0])./norm(p0));
rn = ra*rb/sqrt((ra*sin(psi))^2 + (rb*cos(psi))^2);
r = n(r0,rn)/(sin(theta) * -dndr(r0,rn));
if abs(dot(gV,v0)) == 1 || isnan(theta)
    v1 = v0;
    p1 = p0 + ds.*v1;
else
    a = sign(det([v0' gV']));
    psi = ds * a / r;   
    R = @(psi) [cos(psi) -sin(psi);
        sin(psi) cos(psi)];
    v1 = R(psi)*v0';
    v1 = v1';
    or = p0' + r.*R(a*pi/2)*v0';
    p1 = R(psi)*(p0'-or) + or;
    p1 = p1';
end
alpha = acos(dot(v1,v0));
d = norm(p0-p1);
dserr = abs((d-ds)/ds)*100;
np = n(r0,rn);
end