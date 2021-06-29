% Function for estimating the local radius of curvature of a ray in a
% circular cymetric gradient by Meggott & Meyer-Rochow (1975)

% Input:
% d - distance between n0 and n1
% n0 - outer refractive index
% n1 - inner refractive index
% v - current ray vector
% p - current ray point

% Output:
% r - local radius of curvature for point p

function r = getLocalRadius(n0,n1,d,p,v)
sN = -p';
sN = sN./sqrt(sN(1)^2 + sN(2)^2);
theta = acos(sN(1)*v(1) + sN(2)*v(2));
r = ((n0 + n1)/2) / (sin(theta)*((n1-n0)/d));
end