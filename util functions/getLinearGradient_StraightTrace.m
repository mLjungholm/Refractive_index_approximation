function ns = getLinearGradient_StraightTrace(r0,r1,p0,v0,n0,ns,phase)

[p1,p2,~] = circleIntersect(r0,p0,v0);
d = norm(p2-p1);
% ns = phase/d + n0;

maxsteps = 10^3;
dy = d/maxsteps;
tresh = 0.01;
% n1 = ns;
ex2 = 0;
maxiter = 100;
iter = 1;
while ex2 < tresh && iter < maxiter
    dphase = 0;
    exitvol = 0;
    s = 0;
    p = p1;
    steps = 1;
    while ~exitvol
        k = getK(n0,ns,r0,r1);
        nFunc = @(r) (r-r0)*k + n0;
        r = norm(p);
        if r > r0 && steps > 2
            exitvol = 1;
        elseif steps > maxsteps*2
            exitvol = 1;
        end
        nt = nFunc(r);
        dphase = dphase + dy*(nt-n0);
        p = p + dy.*v0;
        s = s+dy;
        steps = steps + 1;
    end
    acc = abs(phase-dphase);
    if acc > tresh
        if dphase < phase
            ns = ns+0.001;
        else
            ns = ns-0.001;
        end
    else
        ex2 = 1;
    end
    iter = iter + 1;
end

function k = getK(n0,n1,r0,r1)
dn = n1-n0;
dr = r1-r0;
k = dn/dr;
end

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
end