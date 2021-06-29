% Function for the itterative ray trace.

% Input:    Ray positions
%           Initial refractive index
%           Phase shift profile


% Initiate first ray


% Find starting coordinates

% Trace one integration distance
    % Find the local curvature of radius
        % Move to new position and calculate phase shift




function r = findLocalRadius(n0,n1,d,p,v)
% gV = p;
p = p./sqrt(p(1)^2 + p(2)^2);
v = sqrt(v(1)^2 * v(2)^2);
sN = [0 1; -1 0]*p;
theta = acos(sN(1)*v(1) + sN(2)*v(2));
r = ((n0 + n1)/2) / (sin(theta)*((n1-n0)/d));
end

function [p1,v1] = curvePath(p0,v0,r,intDist)
v0 = v0'; p0 = p0';
cv = [0 1; -1 0]*v0;
c = p + r*cv;
psi = intDist/r;
transV = - c;
p1 = ([cos(psi) sin(psi); -sin(psi) cos(psi)]*(p0+transV))-transV;
v1 = [0 1; -1 0]*p1;
p1 = p1'; v1 = v1';
end