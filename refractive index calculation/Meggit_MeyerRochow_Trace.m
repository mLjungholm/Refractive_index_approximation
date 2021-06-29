% Ray tracing based on algorithm by Nilsson, D-E et.al (1983), 
% using an approximation by Meggit & Meyer-Rochow (1975)
% Written by 

function [n,shellR] = Meggit_MeyerRochow_Trace(s,pp,ps,n0,steps)
n = zeros(s.nRays);
n(1) = n0;
shellR = pp(1:end-1);
ds = pp(1)*2/steps;
maxSteps = steps*2;

% First approximation for first ray
[ip0,ip1,intersect] = circleIntersect(pp(1),s.P(1,:),s.V(1,:));
if intersect
    dv = (ip1-ip0);
    d = sqrt(dv(1)^2+dv(2)^2);
    n(2) = n0 + (ps(2)/d);
else
    disp('Error: First aproximation intersect error')
    return
end
% First approximation for the remaining rays
for rayInd = 2:s.nRays
    n(rayInd+1) = firstApproximation(rayInd);
end

% Second approximation for all rays exept last ray
for rayInd = 1:s.nRays-1 % Exclude last ray that needs a special version
    [ip0,~,intersect] = circleIntersect(pp(1),s.P(rayInd,:),s.V(rayInd,:));
    if intersect
        loopRay(ip0,rayInd);
    end
end

    function nr = firstApproximation(rayInd)
        nShells = rayInd; % Number of shells to be taken into acount
        dShells = zeros(nShells,1); % Distance traveled in each shell
        for shellInd = 1:nShells
            [ip0,ip1,~] = circleIntersect(shellR(shellInd),s.P(rayInd,:),s.V(rayInd,:));
            dv = (ip1-ip0);
            d = sqrt(dv(1)^2+dv(2)^2); % distance traveled in shell
            dShells(shellInd) = d;
        end
         % Substract the distance from inner shells.
        for shellInd = 1:nShells-1
            dShells(shellInd) = dShells(shellInd) - dShells(shellInd+1); 
        end
        dPhase = zeros(nShells,1);
        for shellInd = 1:nShells-1
            dPhase(shellInd) = dShells(shellInd)*(n(1+shellInd) - n(1));
        end
        dPhase(rayInd) = ps(rayInd+1) - sum(dPhase(1:rayInd-1));
        nr = dPhase(end)/dShells(end) + n0;
    end

    function loopRay(ip,rayInd)
        maxIter = 100;      % Maximum number of itterations
        threshold = 0.001;  % Threshold value
        iter = 1;           % Current itteration
        accuarcy = inf;     % Current accuarcy
        nRef = 0;           % Current calculated refractive index
        nShells = rayInd;   % Number of shells the ray needs to pass
        rayPath = [];       % Path of ray
        % Run untill threshold is met or max itterations are reached
        while iter < maxIter && accuarcy > threshold
%             pathInShell = zeros(3,1);
            rayPath = zeros(maxSteps,2).*nan;  % Create new path vector
            phasePath = 0;
            phaseSum = 0;
            rayPath(1,:) = s.P(rayInd,:);   % Set origin point
            rayPath(2,:) = ip;              % Set volume intersection point
            rayStep = 2;                    % Current ray step
            v0 = s.V(rayInd,:);             % Current Vector
            p0 = ip;        % Current point
            minR = inf;                     % Current smalest radius
            exitVolume = 0;                 % Exit flag
            cShell = 1;                     % Current shell nr.
            r0 = sqrt(p0(1)^2 + p0(2)^2);   % Current distance from center.
            % Trace the ray through the volume until exit or max steps are
            % reached
            while ~exitVolume && rayStep < maxSteps
                % Find local curvature of radius
                v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
                
                dr = findLocalRadius(n(cShell),n(cShell+1),(shellR(cShell)-shellR(cShell+1)),p0,v0);
                
                [p1,v1] = curvePath(p0,v0,dr,ds); % Move one integration distance
                
                r1 = sqrt(p1(1)^2 + p1(2)^2);  % Current distance from center of volume
                rayStep = rayStep + 1;
                rayPath(rayStep,:) = p1;         % Add step to ray path
                
                if cShell == rayInd
                    phasePath = phasePath + ds;
                else
                    phaseSum = phaseSum + phaseShiftKnownGradient(r0,r1,cShell);
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
            phaseTot = findTotPhase(p0,v0);
            ns  = n(1) + (phaseTot - phaseSum)/phasePath;
            n(1 + rayInd) = (ns+n(1 + rayInd))/2;
            accuarcy = abs(nRef - n(1+rayInd));
            nRef = n(1+rayInd);
            shellR(1+rayInd) = minR;
            iter = iter + 1;
        end
        
        s.PT{rayInd} = rayPath;
%         [k,m] = getSlope(n(rayInd),n(rayInd+1),shellR(rayInd),shellR(rayInd+1));
%         nSlope(rayInd,:)  = [k,m];
        disp(strcat('RayNr: ',num2str(rayInd),', Itterations :',num2str(iter),', Accuarcy: ',num2str(accuarcy),', n = ',num2str(n(1 + rayInd))))
    end


    function r = findLocalRadius(n0,n1,d,p,v)
        sN = -p';
        sN = sN./sqrt(sN(1)^2 + sN(2)^2);
        theta = acos(sN(1)*v(1) + sN(2)*v(2));
        r = ((n0 + n1)/2) / (sin(theta)*((n1-n0)/d));
    end

    function phaseS = phaseShiftKnownGradient(r0,r1,shellInd)
        k = (n(shellInd+1)-n(shellInd))/(shellR(shellInd +1) - shellR(shellInd));
        if r1 > r0
            dr = (r1-r0);
            dn = n(shellInd) - dr*k;
        else
            dr = (r0-r1);
            dn = n(shellInd+1) + dr*k;
        end        
        phaseS = ds*(dn-n(1));
    end

    function [k,m] = getSlope(n0,n1,r0,r1)
        k = (n1-n0)/(r1-r0);
        m = n1-(k*r1);
    end

    function phaseTot = findTotPhase(p0,v0)
        t = -p0(2)/v0(2);
        xip = p0(1) + t*v0(1); % x-axis intersection distance
        % find colsest phaseShift values and interpolate
        if xip > pp(1)
            phaseTot = 0;
            return
        elseif xip < pp(end)
            phaseTot = pp(end);
            return
        end
        [~,indX] = min(abs(pp - xip));
        if xip <= pp(indX)
            pp1 = indX;
            pp2 = indX+1;
        else
            pp1 = indX-1;
            pp2 = indX;
        end
        shellR(rayInd+1) = minR;
        phaseTot = ((ps(pp2)-ps(pp1))/(pp(pp2)-pp(pp1)))*(xip-pp(pp1));
    end
end












