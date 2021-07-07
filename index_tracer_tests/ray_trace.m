function ray_trace(s,ds,n0,n1,rs,nProfile,stepFunc)
steps = round(rs/ds);
for i = 1:s.nRays
    maxSteps = steps*3;
    rayPath = zeros(maxSteps,2).*nan;
    diviation = zeros(maxSteps,1).*nan;
    stepError = zeros(maxSteps,1).*nan;
    rn = zeros(maxSteps,1).*nan;
    phase = 0;
    totalPath = 0;
    
    [p0,~,~] = circleIntersect(rs,s.P(i,:),s.V(i,:));
    v0 = s.V(i,:);
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
        case 'snell'
            while ~exitflag && step < maxSteps
                [p1,v1,alpha,dserr,np,d] = snellDiff_step(ds,p0,v0,n0,n1,rs,nProfile);
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
    end
    s.PT{i} = rayPath;
    s.diviation{i} = diviation;
    s.totalPath(i) = totalPath;
    s.stepError{i} = stepError;
    s.nPath{i} = rn;
    s.phase(i) = phase;
    s.P(i,:) = p1;
    s.V(i,:) = v1;
end
end