% Functions for testing the different solvers for different index profiles

function [p1,v1,alpha,dserr,np,d] = meggit_step(ds,p0,v0,n0,n1,rn,nProfile)
switch nProfile
    case 'linear'
        k = (n0-n1)/rn;
        n = @(r) k*r + n1;
        dndr = @(r) k;
    case 'parabolic'
        k = (n0-n1)/rn^2;
        n = @(r) k*r^2 + n1;
        dndr = @(r) 2*k*r;
    case 'eliptical'
        k = (n1-n0)/rn;
        n = @(r) k*sqrt(rn^2 - r^2) + n0;
        dndr = @(r) -r*k/sqrt(rn^2 - r^2);
    case 'step'
        k = n1-n0;
        n = @(r) n0 + k*(1-sing(r-rn)/2);
        dndr = @(r) 0;
    case 'luneburg'
        R = rn;
        n = @(r) sqrt(2-r^2/R^2);
        dndr = @(r) -r/sqrt(2-r^2/R^2)/R^2;
end

r0 = norm(p0);
gV = -p0;
theta = acos(dot(v0,gV));
r = n(r0)/(sin(theta) * -dndr(r0));
if r == 0
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
% p1 = r.*(R(a*pi/2) - R(a*(pi/2 + psi)))*v0' + p0';
alpha = acos(dot(v1,v0));
% if alpha ~= psi
%     disp('Error: there was an rotation error in the meggit_step')
% end
d = norm(p0-p1);
dserr = abs((d-ds)/ds)*100;
np = n(r0);
% fprintf('The ds step error in "meggit" is: %.4f %% \n',dserr)
end