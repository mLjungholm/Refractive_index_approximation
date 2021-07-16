% Development script for the gradient aproximation algorithm

% Load data or generate new interference data with the
% "createPhaseShift.m" script.

% create a source with variable spacing to match the positions of the
% measured phase values.

testS = Source_2d([peakPos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% lambda = 550*10^-9;
% r = 12*lambda;

n = homogenious_approximation(testS, mPhaseShift, peakPos,n0);

steps = 10^3;
[ng,shellGradient,shellR] = gradient_approximation(testS,peakPos,mPhaseShift,n,steps);

figure(2)
hold on; grid on
plot(peakPos,n,'*')
k = (n0-n1)/r^2;
nFunc = @(r) k.*r.^2 + n1;
plot(truePhasePos,nFunc(truePhasePos),'b')
scatter(shellR,ng,'ro')
legend('Homogenious approx','True gradient','Linear approx')

% err = (abs(nFunc(peakPos) - ng))./nFunc(peakPos);
% figure(3)
% grid on
% plot(peakPos,err)





% Ray tracing based on algorithm by Nilsson, D-E et.al (1983),
% using an approximation by Meggit & Meyer-Rochow (1975)
% Written by Mikael Ljungholm (2021)

function [n,shellGradient,shellR] = gradient_approximation(source,mPhasePos,mPhaseShift,n_homogenious,steps)
% n = zeros(source.nRays+1,1);
% n(1) = n0;
n = n_homogenious;
shellR = mPhasePos;
ds = mPhasePos(1)*2/steps;
maxSteps = steps*2;
shellGradient = zeros(source.nRays,2).*nan;

setGradient(1);

% Second approximation for all rays exept last ray
% for rayInd = 1:s.nRays-1 % Exclude last ray that needs a special version
% for rayInd = 1:s.nRays-1
for rayInd = 1:source.nRays
   stopflag = loopRay(rayInd);
   if stopflag
       break
   end
end


    function stopflag = loopRay(rayInd)
        stopflag = 0;
        % Set up recurcive aproximation
        maxIter = 200;      % Maximum number of itterations
        threshold = 0.001;  % Threshold value
        iter = 1;           % Current itteration
        accuarcy = inf;     % Current accuarcy
        nRef = 0;           % Current calculated refractive index
        nShells = rayInd;   % Number of shells the ray needs to pass
        rayPath = [];       % Path of ray
        % Run untill threshold is met or max itterations are reached
        
%         fprintf('Ray nr: %u \n',rayInd)
        while iter < maxIter && accuarcy > threshold
            rayPath = zeros(maxSteps,2).*nan;  % Create new path vector
            phasePath = 0;
            phaseSum = 0;
            [p0,~,~] = circleIntersect(mPhasePos(1),source.P(rayInd,:),source.V(rayInd,:));
            rayPath(1,:) = p0;              % Set origin point
            v0 = source.V(rayInd,:);             % Current Vector
            rayStep = 1;                    % Current ray step
            minR = inf;                     % Current smalest radius
            cShell = 1;                     % Current shell nr.
            r0 = sqrt(p0(1)^2 + p0(2)^2);   % Current distance from center.
            
            % Trace the ray through the volume until exit or max steps are
            % reached
            shellPath = zeros(maxSteps,2).*nan;
            shellPath(1,:) = p0;
            totalPath = 0;
            phase = 0;
            
            % Initiate trace loop
            exitVolume = 0;     % Exit flag
            while ~exitVolume && rayStep < maxSteps
                % Find local curvature of radius
                v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
                [p1,v1] = intStep(p0,v0,cShell);
                if v1(2) > 0
                    minR = rayPath(1,1);
                    p0 = [shellR(1) 0];
                    v0 = [1,-1];
                    continue
                end
                r1 = sqrt(p1(1)^2 + p1(2)^2);  % Current distance from center of volume
                rayStep = rayStep + 1;
                rayPath(rayStep,:) = p1;         % Add step to ray path
                dt = norm(p1-p0);
%                 totalPath = totalPath + dt;
                totalPath = totalPath + ds;
                phase = phase + phaseShiftKnownGradient(r0,r1,cShell,dt);
                if cShell == rayInd
%                     phasePath = phasePath + dt;
                    phasePath = phasePath + ds;
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
                if r1 <= shellR(cShell+1) && cShell < nShells
                    cShell = cShell + 1;
                elseif r1 > shellR(cShell)
                    cShell = cShell - 1;
                    if cShell < 1
                        exitVolume = 1;
                    end
                end
            end
            
            % Calculate phase shift, back project ray & compare phase shift
            % values
            % [phaseTot,xip]
            phaseTot = findTotPhase(p0,v0);
            n1  = n(1) + (phaseTot - phaseSum)/phasePath;
%             n1 = getLinearGradient_CurvedTrace(shellR(rayInd),shellR(rayInd+1),n(rayInd),n1,(phaseTot - phaseSum),[shellPath(:,1) shellPath(:,2)]);
            n1 = (n(1+rayInd) + n1)/2;
            
            n(1 + rayInd) = n1;
            accuarcy = abs(nRef - n(1+rayInd));
            nRef = n(1+rayInd);
            shellR(1+rayInd) = minR;
            setGradient(rayInd);
            iter = iter + 1;
            if iter == maxIter
                stopflag = 1;
                fprintf('Error: Exceeding maximum itterations! \n Faulty ray ind = %u \n',rayInd)
                return
            end
            
%             figure(1)
%             hold on;
%             axis equal
%             plot(rayPath(:,1),rayPath(:,2))
%             plotCircle(shellR(1),pi,1)
%             plotLine([0 0],[shellR(1) 0],'k')
%             plotLine([0 shellR(1)],[0 -shellR(1)],'k')
            
%             fprintf('Itteration : %u \n',iter-1);
%             fprintf('n = %.4f \n',n1)
            
        end
        if n(1+rayInd) < 0
            fprintf('Error: approximated refractive index was negative! \n Faulty rayInd = %u \n',rayInd)
            n(1+rayInd) = n_homogenious(1+rayInd);
            stopflag = 1;
            return
        end
        source.P(rayInd,:) = p1;
        source.V(rayInd,:) = v1;
        source.path{rayInd} = rayPath;
        source.phase(rayInd) = phase;
        source.totalPath(rayInd) = totalPath;
    end

    function [p1,v1] = intStep(p0,v0,shellInd)
        r0 = norm(p0);
        gV = -p0./r0;
        gV = gV./norm(gV);
        theta = acos(dot(v0,gV)./norm(v0));
        nr = shellGradient(shellInd,1)*(r0-shellR(shellInd)) + shellGradient(shellInd,2);
        r = nr/(sin(theta) * -shellGradient(shellInd,1));
        if abs(dot(gV,v0)) == 1 || isnan(theta)
            v1 = v0;
            p1 = p0 + ds.*v1;
        else
            a = sign(det([v0' gV']));
            psi = ds * a / r;
            R = @(psi) [cos(psi) -sin(psi);
                sin(psi) cos(psi)];
            v1 = R(psi)*v0';
            v1 = v1';
            or = p0' + r.*R(a*pi/2)*v0';
            p1 = R(psi)*(p0'-or) + or;
            p1 = p1';
        end
    end

    function phaseS = phaseShiftKnownGradient(r0,r1,shellInd,dt)
        nr = shellGradient(shellInd,1)*((r0+r1)/2-shellR(shellInd)) + shellGradient(shellInd,2);
        phaseS = dt*(nr-n(1));
    end

    function setGradient(shellInd)
        if isnan(shellGradient(1,1))
            for si = 1:size(shellGradient,1)
                [k,m] = getGradient(si);
                shellGradient(si,:) = [k,m];
            end
        elseif shellInd == size(shellGradient,1)
            [k,m] = getGradient(shellInd);
            shellGradient(shellInd,:) = [k,m];
        else
            [k,m] = getGradient(shellInd);
            shellGradient(shellInd,:) = [k,m];
            [k,m] = getGradient(shellInd+1);
            shellGradient(shellInd+1,:) = [k,m];
        end
        function [k,m] = getGradient(si)
            na = n(si);
            nb = n(si+1);
            ra = shellR(si);
            rb = shellR(si+1);
            dn = na-nb;
            dr = ra-rb;
            k = dn/dr;
            m = na;
        end
    end

    function phaseTot = findTotPhase(p0,v0)
        t = -p0(2)/v0(2);
        xip = p0(1) + t*v0(1); % x-axis intersection distance
        % find colsest phaseShift values and interpolate
        pind = find((mPhasePos-xip)>0,1,'last');
        if isempty(pind)
            pind = 1;
        elseif pind == length(mPhasePos)
            pind = length(mPhasePos) - 1;
        end
        p0 = mPhasePos(pind);
        p1 = mPhasePos(pind+1);
        shift0 = mPhaseShift(pind);
        shift1 = mPhaseShift(pind+1);       
        phaseTot = (shift1-shift0)/(p1-p0)*(xip-p0) + shift0;
        if phaseTot < 0
            phaseTot = 0;
        end
    end

end














