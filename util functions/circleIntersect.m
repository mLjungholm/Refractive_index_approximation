function [p1,p2,intersect] = circleIntersect(r,p0,v)
L = -p0;
tca = v(1)*L(1) + v(2)*L(2);
if tca < 0
    intersect = 0;
    p1 = nan;
    p2 = nan;
    return
end
d = sqrt(L(1)^2 + L(2)^2 - tca^2);
if d < 0 || d > r
    intersect = 0;
    p1 = nan;
    p2 = nan;
    return
end
thc = sqrt(r^2-d^2);
t0 = tca - thc;
t1 = tca + thc;
p1 = p0 + t0*v;
p2 = p0 + t1*v;
intersect = 1;
end