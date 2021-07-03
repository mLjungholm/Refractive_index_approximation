function [p1,v1,alpha,dserr,np,d] = snellDiff_step(ds,p0,v0,n0,n1,rn,nProfile)
switch nProfile
    case 'linear'
        k = (n0-n1)/rn;
        n = @(r) k*r + n1;
%         dndr = @(r) k;
    case 'parabolic'
        k = (n0-n1)/rn^2;
        n = @(r) k*r^2 + n1;
%         dndr = @(r) 2*k*r;
    case 'eliptical'
        % Ingore for now
    case 'step'
        k = n1-n0;
        n = @(r) n0 + k*(1-sing(r-rn)/2);
%         dndr = @(r) 0;
end

r0 = norm(p0);
gV = p0;
gV = gV./norm(gV);
% theta = acos(dot(v0,gV));
ra = r0 - ds/2.*gV;
rb = r0 + ds/2.*gV;
na = n(ra);
nb = n(rb);

dv = dot(v0,gV);
if dv < 0
    gV = -gV;
end

theta = acos(dv);


% if dv < 0
% theta = theta -pi/2;
% end
% a = na/nb*sin(theta);
if theta == 0
    v1 = v0;
    p1 = p0 + ds.*v0;
elseif theta == pi/2
    v1 = v0;
    p1 = p0 + ds.*v0;
elseif a > 1
    v1 = -v0;
    p1 = p0 + ds.*v1;
else
    psi = theta - asin(a);
    psi = b*psi;
    R = [cos(psi) -sin(psi);
            sin(psi) cos(psi)];
        v1 = (R*v0')';
        p1 = p0 + ds.*v1;
end
    
% psi = asin(na

% c = dot(v0,gV);
% c = -c*sign(c);
% r = na/nb;
% v1 = r*v0 + (r*c - sqrt(1- r^2*(1-c^2)))*gV.*ds;
% v1 = v1./norm(v1);
% p1 = p0 + ds.*v1;

if ~isreal(v1)  
    disp('Error: total internal reflection occured in snellDiff_step')
end
alpha = acos(dot(v0,v1));
d = norm(p0-p1);
dserr = (d-ds)/ds*100;
np = n(r0);
% fprintf('The ds step error in "meggit" is: %.4f %% \n',dserr)
end