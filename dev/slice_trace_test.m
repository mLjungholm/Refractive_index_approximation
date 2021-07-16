function slice_trace_test(s,steps,R,n0,k,rayInd)
ds = R/steps;
maxSteps = steps*2;
rayPath = zeros(maxSteps,2).*nan;  % Create new path vector
phasePath = 0;
phaseSum = 0;
[p0,~,~] = circleIntersect(R,s.P(rayInd,:),s.V(rayInd,:));
rayPath(1,:) = p0;              % Set origin point
v0 = s.V(rayInd,:);             % Current Vector
rayStep = 2;                    % Current ray step
minR = inf;                     % Current smalest radius
cShell = 1;                     % Current shell nr.
r0 = sqrt(p0(1)^2 + p0(2)^2);   % Current distance from center.
% Trace the ray through the volume until exit or max steps are reached
shellPath = zeros(maxSteps,2).*nan;
shellPath(1,:) = p0;
totalPath = 0;
phase = 0;

% Initiate trace loop
exitVolume = 0;     % Exit flag
while ~exitVolume && rayStep < maxSteps
    v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
    [p1,v1] = intStep(p0,v0,cShell);
    if v1(2) > 0 % Catch to high refractive index. I.e if the ray loops around
        minR = rayPath(1,1);
        p0 = [shellR(1) 0];
        v0 = [1,-1];
        continue
    end
    r1 = sqrt(p1(1)^2 + p1(2)^2);  % Current distance from center of volume
    rayStep = rayStep + 1;
    rayPath(rayStep,:) = p1;         % Add step to ray path
    dt = norm(p1-p0);
    totalPath = totalPath + dt;
    phase = phase + phaseShiftKnownGradient(r0,r1,cShell,dt);
    if cShell == rayInd
        phasePath = phasePath + dt;
        shellPath(rayStep,:) = p1;
    else
        phaseSum = phaseSum + phaseShiftKnownGradient(r0,r1,cShell,dt);
    end
    v0 = v1;                         % Set new vector
    p0 = p1;                         % Set new curent point
    r0 = r1;
    if r1 < minR
        minR = r1;   % New min distance to center of colume
    end
    if r1 > R
        exitVolume = 1;
    end
end
s.P(rayInd,:) = p1;
s.V(rayInd,:) = v1;
s.path{rayInd} = rayPath;
s.phase(rayInd) = phase;
s.totalPath(rayInd) = totalPath;

    function [p1,v1] = intStep(p0,v0,~)
        r0 = norm(p0);
        gV = -p0./r0;
        gV = gV./norm(gV);
        theta = acos(dot(v0,gV)./norm(v0));
        nr = k*(r0-R) + n0;
        %         nr = shellGradient(shellInd,1)*(r0-shellR(shellInd)) + shellGradient(shellInd,2);
        r = nr/(sin(theta) * -k);
        if abs(dot(gV,v0)) == 1 || isnan(theta)
            v1 = v0;
            p1 = p0 + ds.*v1;
        else
            a = sign(det([v0' gV']));
            psi = ds * a / r;
            Rot = @(psi) [cos(psi) -sin(psi);
                sin(psi) cos(psi)];
            v1 = Rot(psi)*v0';
            v1 = v1';
            or = p0' + r.*Rot(a*pi/2)*v0';
            p1 = Rot(psi)*(p0'-or) + or;
            p1 = p1';
        end
    end

    function phaseS = phaseShiftKnownGradient(r0,r1,~,dt)
        nr = k*((r0+r1)/2-R) + n0;
        %         nr = shellGradient(shellInd,1)*((r0+r1)/2-shellR(shellInd)) + shellGradient(shellInd,2);
        phaseS = dt*(nr-n0);
    end
end