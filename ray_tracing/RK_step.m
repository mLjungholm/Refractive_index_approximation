
function [p1,v1,alpha,dserr,np,d] = RK_step(ds,p0,v0,n0,n1,rn,nProfile)
switch nProfile
    case 'parabolic'
        k = (n0-n1)/rn^2;
        n = @(r) k*(r(1)^2 + r(2)^2) + n1;
        D = @(r) 2*k*(k*(r(1)^2 + r(2)^2) + n1).*[r(1) r(2)];
    case 'linear'
        k = (n0-n1)/rn;
        n = @(r) k*sqrt(r(1)^2 + r(2)^2) + n1;
        dndr = @(r) k/sqrt(r(1)^2 + r(2)^2).*[r(1) r(2)];
        D = @(r) k*(k + n1/sqrt(r(1)^2 + r(2)^2)).*[r(1) r(2)];
    case 'step'
        k = n1-n0;
        n = @(r) n0 + k*(1-sign(sqrt(r(1)^2 + r(2)^2)-rn)/2);
        D = @(r) 0;
    case 'eliptical'
        k = (n1-n0)/rn;
        n = @(r) k*sqrt(rn^2 - r(1)^2 - r(2)^2) + n0;
%         dndr = @(r) -k/sqrt(rn^2 - r(1)^2 - r(2)^2).*[r(1) r(2)];
        D = @(r) -(k*sqrt(rn^2 - r(1)^2 - r(2)^2) + n0).*k./sqrt(rn^2 - r(1)^2 - r(2)^2).*[r(1) r(2)];
    case 'luneburg'
        R = rn;
        n = @(r) sqrt(2-(r(1)^2 + r(2)^2)/R^2);
        dndr = @(r) -1/sqrt(2-(x^2 + y^2)/R^2)/R^2.*[x y];
        D = @(r) -1/R^2.*[r(1) r(2)];
%         n = @(r) sqrt(2-r.^2/R^2);
%         D = @(r) -1/R^2.*r;
end
T = v0.*n(p0);
A = real(ds.*D(p0));
if abs(A(1)) > 1
    v1 = v0;
    p1 = p0 + ds.*v1;
else
    B = real(ds.*D(p0 + ds/2.*T + ds/8.*A));
    C = real(ds.*D(p0 + ds.*T + ds/2.*B));
    p1 = p0 + ds.*(T + 1/6.*(A + 2.*B));
%     if ~isreal(p1)
%         testflag = 1;
%     end
    T1 = T + 1/6.*(A + 4.*B + C);
%     v1 = v1./sqrt(v1(1)^2 + v1(2)^2);
%     p1 = p0 + ds.*v1;
v1 = T1./n(p1);
end
alpha = acos(dot(v1,v0)./(norm(v1)*norm(v0)));
d = norm(p0-p1);
dserr = abs((d-ds)/ds)*100;
np = n(p0);
% fprintf('The ds step error in "Runge-Kutta" is: %.4f %% \n',dserr)
end