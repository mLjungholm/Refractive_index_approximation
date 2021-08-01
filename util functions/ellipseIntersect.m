function [p1,p2,intersect] = ellipseIntersect(ra,rb,p0,v)
% 
% x1 = p0(1);
% x2 = p0(1) + v(1);
% y1 = p0(2);
% y2 = p0(2) + v(2);

% A = (x2-x1)^2/ra^2 + (y2-y1)^2/rb^2;
% B = 2*x1*(x2-x1)/ra^2 + 2*y1*(y2-y1)/rb^2;
% C = x1^2/ra^2 + y1^2/rb^2 - 1;

A = v(1)^2/ra^2 + v(2)^2/rb^2;
B = 2*p0(1)*v(1)/ra^2 + 2*p0(2)*v(2)/rb^2;
C = p0(1)^2/ra^2 + p0(2)^2/rb^2 - 1;

D = B^2-4*A*C;
% D2 = B2^2-4*A2*C2;

if D < 0
    intersect = 0;
    p1 = nan;
    p2 = nan;
    return
elseif D == 0
    intersect = 1;
    t = -B/2/A;
    p1 = p0 + v.*t;
    p2 = p1;
    return
end

t1 = (-B - sqrt(D))/2/A;
t2 = (-B + sqrt(D))/2/A;
p1 = p0 + v.*t1;
p2 = p0 + v.*t2;
% p1 = [x1 + (x2-x1)*t1, y1 + (y2-y1)*t1];
% p2 = [x1 + (x2-x1)*t2, y1 + (y2-y1)*t2];
intersect = 1;
end