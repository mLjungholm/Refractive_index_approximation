% Simple function for finding the refractive index gradient in one shell
% that results in the same phase shift as the average refractive index "ns"
function n1 = getLinearGradient_CurvedTrace(r0,r1,n0,ns,phase,rayPath)
rayPath = rayPath(~isnan(rayPath(:,1)),:);
if phase <= 0
    n1 = n0;
    return
end
tresh = 0.01;
ex2 = 0;
maxiter = 100;
iter = 1;
n1 = ns;
while ex2 < tresh && iter < maxiter
    dphase = 0;
    dn = n1-n0;
    dr = r1-r0;
    k = dn/dr;
    nFunc = @(r) (r-r0)*k + n0;
    for rayInd = 1:length(rayPath(:,1))-1
        r = sqrt(rayPath(rayInd,1)^2 + rayPath(rayInd,2)^2);
        nt = nFunc(r);
        ds = norm(rayPath(rayInd+1)-rayPath(rayInd));
        dphase = dphase + ds*(nt-n0);
    end
    acc = abs(phase-dphase);
    if acc > tresh
        if dphase < phase
            n1 = n1+0.001;
        else
            n1 = n1-0.001;
        end
    else
        ex2 = 1;
    end
    iter = iter + 1;
end
end