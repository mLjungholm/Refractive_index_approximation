% Function for moving a point "p0" a distance "ds" along the arc of the
% circle with radius "r" with a radius vector perpendicular to the vector
% "v0"

% Input:
% p0 - starting point
% v0 - starting vector
% r - radius of curvature
% ds - arclength

% Output:
% p1 - end point
% v1 - end vector

function [p1,v1] = curvePath(p0,v0,r,ds)
v0 = v0'; p0 = p0';
detCV = -p0(1)*v0(2) + p0(2)*v0(1);
if detCV > 0
    rotM = [0 1; -1 0];
else
    rotM = [0 -1; 1 0];
end
cv = rotM*v0;
c = p0 + r*cv;
if detCV > 0
    psi = ds/r;
else
    psi = -ds/r;
end
transV = - c;
p1 = p0 + transV;
rotMat =  [cos(psi) sin(psi); -sin(psi) cos(psi)];
p1 = rotMat*p1;
p1 = p1 - transV;
v1 = rotMat*v0;
p1 = p1'; v1 = v1';
end