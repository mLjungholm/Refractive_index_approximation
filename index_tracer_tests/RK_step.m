
function [p1,v1,alpha,dserr,np,d] = RK_step(ds,p0,v0,n0,n1,rn,nProfile)
switch nProfile
    case 'parabolic'
        k = (n0-n1)/rn^2;
        n = @(r) k*(r(1)^2 + r(2)^2) + n1;
        D = @(r) 2*k*(k*(r(1)^2 + r(2)^2) + n1).*[r(1) r(2)];
    case 'linear'
        k = (n0-n1)/rn;
        n = @(r) k*sqrt(r(1)^2 + r(2)^2) + n1;
        D = @(r) k*(k + n1/sqrt(r(1)^2 + r(2)^2)).*[r(1) r(2)];
    case 'step'
        k = n1-n0;
        n = @(r) n0 + k*(1-sing(r-rn)/2);
        D = @(r) 0;
    case 'eliptical'
        %ignore for now
end

A = ds.*D(p0);
B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
C = ds.*D(p0 + ds.*v0 + ds/2.*B);
p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
v1 = v0 + 1/6.*(A + 4.*B + C);
v1 = v1./sqrt(v1(1)^2 + v1(2)^2);

alpha = acos(dot(v1,v0));
d = norm(p0-p1);
dserr = (d-ds)/ds*100;
np = n(p0);
% fprintf('The ds step error in "Runge-Kutta" is: %.4f %% \n',dserr)
end