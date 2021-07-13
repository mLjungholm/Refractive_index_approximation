% Ray tracing based on algorithm by Nilsson, D-E et.al (1983), 
% using an approximation by Meggit & Meyer-Rochow (1975)
% Written by Mikael Ljungholm (2021)

function [n,shellGradient,shellR] = findRefractiveIndex(s,pp,ps,n0,steps)
n = zeros(s.nRays+1,1);
n(1) = n0;
shellR = pp;
ds = pp(1)*2/steps;
maxSteps = steps*2;
shellGradient = zeros(s.nRays,2).*nan;

% First approximation for first ray
[ip0,ip1,intersect] = circleIntersect(pp(1),s.P(1,:),s.V(1,:));
if intersect
    dv = (ip1-ip0);
    d = norm(dv);
    n(2) = n0 + (ps(2)/d);
else
    disp('Error: Ray did not intersect volume')
    return
end

% First approximation for the remaining rays
for rayInd = 2:s.nRays
    n(rayInd+1) = firstApproximation(rayInd);
end
nf = n;
setGradient(1);
s.resetSource();

% Second approximation for all rays exept last ray
% for rayInd = 1:s.nRays-1 % Exclude last ray that needs a special version
% for rayInd = 1:s.nRays-1
for rayInd = 1:1
    loopRay(rayInd);
end

    % Calculates the first approximation of the refractive index ine each
    % shell. This assumes that the rays travel straight threw the slice.
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
        dPhase = zeros(nShells-1,1);
        for shellInd = 1:nShells-1
            dPhase(shellInd) = dShells(shellInd)*(n(1+shellInd) - n(1));
        end
        dPhase(rayInd) = ps(rayInd+1) - sum(dPhase);
        nr = dPhase(end)/dShells(end) + n0;
    end

    function loopRay(rayInd)
        % Set up recurcive aproximation
        maxIter = 200;      % Maximum number of itterations
        threshold = 0.001;  % Threshold value
        iter = 1;           % Current itteration
        accuarcy = inf;     % Current accuarcy
        nRef = 0;           % Current calculated refractive index
        nShells = rayInd;   % Number of shells the ray needs to pass
        rayPath = [];       % Path of ray
        % Run untill threshold is met or max itterations are reached
                
        disp(strcat('Ray nr :',num2str(rayInd)))
        while iter < maxIter && accuarcy > threshold
%           % Initiate ray
            rayPath = zeros(maxSteps,2).*nan;  % Create new path vector
            phasePath = 0;
            phaseSum = 0;
            [p0,~,~] = circleIntersect(pp(1),s.P(rayInd,:),s.V(rayInd,:));
            rayPath(1,:) = p0;              % Set origin point
            v0 = s.V(rayInd,:);             % Current Vector
            rayStep = 2;                    % Current ray step
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
%                 dr = findLocalRadius(n(cShell),n(cShell+1),(shellR(cShell)-shellR(cShell+1)),p0,v0);
%                 if dr < 0
%                     errflag = 1;
%                 end
%                 [p1,v1] = curvePath(p0,v0,dr,ds); % Move one integration distance
                [p1,v1] = intStep(p0,v0,cShell);
                if v1(2) > 0
%                     errflag = 1;
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
            [phaseTot,~] = findTotPhase(p0,v0);
%             if phaseTot == 0
%                 n(rayInd + 1) = n(rayInd + 1)*0.999;  
%                 continue
%             elseif phaseTot <= phaseSum
%                 n(rayInd + 1) = n(rayInd + 1)*0.999;
%             end
            n1  = n(1) + (phaseTot - phaseSum)/phasePath;
%             n1 = getLinearGradient_CurvedTrace(shellR(rayInd),shellR(rayInd+1),n(rayInd),ns,(phaseTot - phaseSum),[shellPath(:,1) shellPath(:,2)],ds);
            n1 = (n(1+rayInd) + n1)/2;
            
            n(1 + rayInd) = n1;
            accuarcy = abs(nRef - n(1+rayInd));
            nRef = n(1+rayInd);
            shellR(1+rayInd) = minR;
            setGradient(rayInd);
            iter = iter + 1;
            
            figure(1)
            hold on;
            axis equal
            plot(rayPath(:,1),rayPath(:,2))
            plotCircle(shellR(1),pi,1)
%             plotCircle(shellR(rayInd+1),pi,1)
            plotLine([0 0],[shellR(1) 0],'k',1)
            plotLine([0 shellR(1)],[0 -shellR(1)],'k',1)
            
            fprintf('Itteration : %u \n',iter-1);
            
        end
        figure(1)
%         hold on; 
%         axis equal
% %         yyaxis left
%         plot(rayPath(:,1),rayPath(:,2),'b')
%         plotCircle(shellR(1),pi,1)
        plotCircle(shellR(rayInd+1),pi,1)
%         plotLine([0 0],[shellR(1) 0],'k',1)
%         plotLine([0 shellR(1)],[0 -shellR(1)],'k',1)
        
        s.P(rayInd,:) = p1;
        s.V(rayInd,:) = v1;
        s.PT{rayInd} = rayPath;
        s.phase(rayInd) = phase;
        s.totalPath(rayInd) = totalPath;
    end


%     function r = findLocalRadius(n0,n1,d,p,v)
%         sN = -p';
%         sN = sN./sqrt(sN(1)^2 + sN(2)^2);
%         theta = acos(sN(1)*v(1) + sN(2)*v(2));
%         r = ((n0 + n1)/2) / (sin(theta)*((n1-n0)/d));
%     end

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

    function [phaseTot,xip] = findTotPhase(p0,v0)
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
        phaseTot = ((ps(pp2)-ps(pp1))/(pp(pp2)-pp(pp1)))*(xip-pp(pp1));
    end
end













