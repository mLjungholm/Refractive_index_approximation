function ray_trace_single_ray(s,steps,n0,n1,rs,nProfile,stepFunc,rayNr)
% steps = round(rs/ds);
ds = rs/steps;

maxSteps = steps*3;
rayPath = zeros(maxSteps,2).*nan;
diviation = zeros(maxSteps,1).*nan;
stepError = zeros(maxSteps,1).*nan;
rn = zeros(maxSteps,1).*nan;
phase = 0;
totalPath = 0;

[p0,p02,intersect] = circleIntersect(rs,s.P(rayNr,:),s.V(rayNr,:));
if ~intersect || isequal(p0, p02)
    s.path{rayNr} = [0 0];
    s.diviation{rayNr} = nan;
    s.totalPath(rayNr) = nan;
    s.stepError{rayNr} = nan;
    s.phase(rayNr) = nan;
    return
end
v0 = s.V(rayNr,:);
step = 1;
exitflag = 0;
switch stepFunc
    case 'meggit'
        while ~exitflag && step < maxSteps
            [p1,v1,alpha,dserr,np,d] = meggit_step(ds,p0,v0,n0,n1,rs,nProfile);
            r1 = norm(p1);
            if r1 > rs
                exitflag = 1;
            end
            rayPath(step+1,:) = p1;
            phase = phase + norm(p1-p0)*(np-n0);
            totalPath = totalPath + d;
            diviation(step) = alpha;
            stepError(step) = dserr;
            rn(step) = np;
            p0 = p1;
            v0 = v1;
            step = step + 1;
        end
    case 'megg_diff'
        while ~exitflag && step < maxSteps
            [p1,v1,alpha,dserr,np,d] = meggit_diff_step(ds,p0,v0,n0,n1,rs,nProfile);
            r1 = norm(p1);
            if r1 > rs
                exitflag = 1;
            end
            rayPath(step+1,:) = p1;
            phase = phase + norm(p1-p0)*(np-n0);
            totalPath = totalPath + d;
            diviation(step) = alpha;
            stepError(step) = dserr;
            rn(step) = np;
            p0 = p1;
            v0 = v1;
            step = step + 1;
        end
    case 'rk'
        while ~exitflag && step < maxSteps
            [p1,v1,alpha,dserr,np,d] = RK_step(ds,p0,v0,n0,n1,rs,nProfile);
            r1 = norm(p1);
            if r1 > rs
                exitflag = 1;
            end
            rayPath(step+1,:) = p1;
            phase = phase + norm(p1-p0)*(np-n0);
            totalPath = totalPath + d;
            diviation(step) = alpha;
            stepError(step) = dserr;
            rn(step) = np;
            p0 = p1;
            v0 = v1;
            step = step + 1;
        end
    case 'meggitStepGradient'
        shell = 1;
        while ~exitflag && step < maxSteps
            n1 = nProfile(shell,2);
            [p1,v1,alpha,dserr,np,d] = meggit_step(ds,p0,v0,nProfile(shell,1),n1,nProfile(shell,3),'linear');
            r1 = norm(p1);
            if r1 > rs
                exitflag = 1;
            elseif r1 < nProfile(shell+1,3)
                shell = shell + 1;
            elseif r1 > nProfile(shell,3)
                shell = shell - 1;
            end
            rayPath(step+1,:) = p1;
            phase = phase + norm(p1-p0)*(np-n0);
            totalPath = totalPath + d;
            diviation(step) = alpha;
            stepError(step) = dserr;
            rn(step) = np;
            p0 = p1;
            v0 = v1;
            step = step + 1;
        end
end
s.path{rayNr} = rayPath;
s.diviation{rayNr} = diviation;
s.totalPath(rayNr) = totalPath;
s.stepError{rayNr} = stepError;
s.phase(rayNr) = phase;
s.P(rayNr,:) = p1;
s.V(rayNr,:) = v1;

end