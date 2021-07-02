function meggitMeyer_rayTrace_known_gradient(s,r,steps,nRange)
nRays = s.nRays;
ds = r*2/steps;
maxSteps = steps*2;
n0 = nRange(1);
n1 = nRange(2);
k = (n0-n1);
nFunc = @(x) k.*x.^2 + n1; 
arrayfun(@traceRay,1:s.nRays);

% function itterativeTrace_singleGradient(s,r,steps)
    function traceRay(rayInd)

% for rayInd = 1:s.nRays
    [ip0,ip1,intersect] = circleIntersect(r,s.P(rayInd,:),s.V(rayInd,:));
    if intersect
        if ip0 == ip1
            rayPath = zeros(2,2);
            rayPath(1,:) = s.P(rayInd,:);   % Set origin point
            yStop = -1;
            endP = [ip0(1),yStop];
            rayPath(2,:) = endP;
            s.P(rayInd,:) = endP;
            s.PT{rayInd} = rayPath;
            s.phase(rayInd) = 0;
            return
%             s.phase(rayInd) = norm(endP-s.P0(rayInd,:))*s.nStart;
%             continue
        end
        loopRay(ip0,rayInd);
    else
        rayPath = zeros(2,2);
        yStop = -1;
        endP = [s.P(rayInd,1),yStop];
        rayPath(1,:) = s.P(rayInd,:);
        rayPath(2,:) = endP;
        s.P(rayInd,:) = endP;
        s.PT{rayInd} = rayPath;
        s.phase(rayInd) = 0;
    end
end

    function loopRay(ip,rayInd)
        rayPath = [];
        exitVol = 0;
        rayPath = zeros(maxSteps,2);  % Create new path vector
        vPath = rayPath.*nan;
        nPath = rayPath.*nan;
        rayPath(:,1) = nan; rayPath(:,2) = nan;
        rayPath(1,:) = s.P(rayInd,:);   % Set origin point
        rayPath(2,:) = ip;              % Set volume intersection point
        %         phaseSum = norm(ip-s.P0(rayInd,:))*s.nStart;
        phaseSum = 0;
        rayStep = 2;                    % Current ray step
        v0 = s.V(rayInd,:);             % Current Vector
        v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
        p0 = rayPath(rayStep,:);        % Current point
        r0 = sqrt(p0(1)^2 + p0(2)^2);   % Current distance from center.
        vPath(1,:) = v0;
        nPath(1) = n0;
        % Run untill threshold is met or max itterations are reached
        while ~exitVol && rayStep < maxSteps
            %             sN = -p0';
            %             sN = sN./sqrt(sN(1)^2 + sN(2)^2);
            %             theta = acos(sN(1)*v0(1) + sN(2)*v0(2));
            %             rn = ((n0 + n1)/2) / (sin(theta)*(n1-n0));
            vPath(rayStep,:) = v0;
            nPath(rayStep) = nFunc(r0);
            rn = getLocalRadiusKnownGradient(p0,v0,r0);
            [p1,v1] = curvePath(p0,v0,rn,ds); % Move one integration distance           
            r1 = sqrt(p1(1)^2 + p1(2)^2);  % Current distance from center of volume
%             if ~isreal(r1)
%                 testflag = 1;
%             end           
            if r1 > r && rayStep > 2 
                yStop = -1;
                tStop = (yStop-p1(2))/v1(2);
                xStop = p1(1) + tStop*v1(1);
                p0 = p1;
                rayPath(rayStep+1,:) = p1;
                p1 = [xStop,yStop];
%                 phaseSum = phaseSum + norm(p1-p0));
                exitVol = 1;
            end
            
            rayStep = rayStep + 1;
            rayPath(rayStep,:) = p1;         % Add step to ray path
%             phaseSum = phaseSum + ds*abs(nGradient(r0)-s.nStart);
            phaseSum = phaseSum + ds*abs(nFunc(r0)-s.nStart);
%             r0 = r1;
            v0 = v1;                         % Set new vector
            p0 = p1;    % Set new curent point
            r0 = r1;
        end
        s.PT{rayInd} = rayPath;
        s.nPath{rayInd} = nPath;
        s.vPath{rayInd} = vPath;
        s.phase(rayInd) = phaseSum;
        s.P(rayInd,:) = p0;
        s.V(rayInd,:) = v0;
%     end
end



function r = getLocalRadiusKnownGradient(p,v,r0)
% nGradient = @(r) real(sqrt(1 - r.^2)).*(1.3-1) + 1;
sN = -p';
sN = sN./sqrt(sN(1)^2 + sN(2)^2);
% rs = sqrt(p(1)^2 + p(2)^2);
% ns = nFunc(rs);
dn = nFunc(r0);
% na = nFunc(r0 + ds/2);
% nb = nFunc(r0 - ds/2);
% if na == 1
%     D = norm(p-sN.*ds) - 1;
% else
%     D = 2*ds;
% end
% D = ds;
% D = sqrt(r0^2 + r1^2);
theta = acos(sN(1)*v(1) + sN(2)*v(2));
% r = ((na + nb)/2) / (sin(theta)*((nb-na)/(D)));
r = (dn) / (sin(theta)*(-2*k*r0)*2);
test = 1;
end
end