function ray_trace_elliptical(s,steps,n0,n1,ra,rb,nProfile,stepFunc)
% steps = round(rs/ds);
ds = max([ra,rb])/steps;
for i = 1:s.nRays
    maxSteps = steps*3;
    rayPath = zeros(maxSteps,2).*nan;
    diviation = zeros(maxSteps,1).*nan;
    stepError = zeros(maxSteps,1).*nan;
    rn = zeros(maxSteps,1).*nan;
    phase = 0;
    totalPath = 0;
    
    [p0,p02,intersect] = ellipseIntersect(ra,rb,s.P(i,:),s.V(i,:));
    if ~intersect || isequal(p0, p02)
        s.path{i} = [nan nan];
        s.diviation{i} = nan;
        s.totalPath(i) = nan;
        s.stepError{i} = nan;
        %     s.nPath{i} = rn;
        s.phase(i) = nan;
%         s.P(i,:) = p1;
%         s.V(i,:) = v1;
        continue
    end
    rayPath(1,:) = p0;
    v0 = s.V(i,:);
    step = 1;
    exitflag = 0;
    switch stepFunc
        case 'meggit'
            while ~exitflag && step < maxSteps
                [p1,v1,alpha,dserr,np,d] = meggit_step_elliptical(ds,p0,v0,n0,n1,ra,rb,nProfile);
%                 r1 = norm(p1);
                rs = p1(1)^2/ra^2 + p1(2)^2/rb^2;
                if rs > 1
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
%%%%%%%%%%%%% These have not been modified to work with elliptical shapes
%         case 'megg_diff'
%             while ~exitflag && step < maxSteps
%                 [p1,v1,alpha,dserr,np,d] = meggit_diff_step(ds,p0,v0,n0,n1,rs,nProfile);
%                 r1 = norm(p1);
%                 if r1 > rs
%                     exitflag = 1;
%                 end
%                 rayPath(step+1,:) = p1;
%                 phase = phase + norm(p1-p0)*(np-n0);
%                 totalPath = totalPath + d;
%                 diviation(step) = alpha;
%                 stepError(step) = dserr;
%                 rn(step) = np;
%                 p0 = p1;
%                 v0 = v1;
%                 step = step + 1;
%             end
%         case 'rk'
%             while ~exitflag && step < maxSteps
%                 [p1,v1,alpha,dserr,np,d] = RK_step(ds,p0,v0,n0,n1,rs,nProfile);
%                 r1 = norm(p1);
%                 if r1 > rs
%                     exitflag = 1;
%                 end
%                 rayPath(step+1,:) = p1;
%                 phase = phase + norm(p1-p0)*(np-n0);
%                 totalPath = totalPath + d;
%                 diviation(step) = alpha;
%                 stepError(step) = dserr;
%                 rn(step) = np;
%                 p0 = p1;
%                 v0 = v1;
%                 step = step + 1;
%             end
%         case 'snell'
%                 N = p0./norm(p0);
%                 [v1, ~] = Snell(v0, N, n0, n1);
%                 [~,p1,~] = circleIntersect(rs,p0,v1);
%                 rayPath(1,:) = p0;                
%                 rayPath(2,:) = p1;
%                 N = p1./norm(p1);
%                 [v1, ~] = Snell(v1, N, n1, n0);
%                 totalPath = norm(p0-p1);
%                 phase = totalPath*(n1-n0);
%                 diviation(step) = 0;
%                 stepError(step) = 0;
% %                 rn(step) = np;
% %                 p0 = p1;
% %                 v0 = v1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    s.path{i} = rayPath;
    s.diviation{i} = diviation;
    s.totalPath(i) = totalPath;
    s.stepError{i} = stepError;
%     s.nPath{i} = rn;
    s.phase(i) = phase;
    s.P(i,:) = p1;
    s.V(i,:) = v1;
end
end