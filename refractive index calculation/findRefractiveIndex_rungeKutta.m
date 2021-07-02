% Ray tracing based on algorithm by Nilsson, D-E et.al (1983), 
% using an approximation by Meggit & Meyer-Rochow (1975)
% Written by Mikael Ljungholm (2021)

function [n,nf,shellR] = findRefractiveIndex_rungeKutta(s,phasePos,phaseVal,n0,steps)
n = zeros(s.nRays,1);
n(1) = n0;
shellR = phasePos;
ds = phasePos(1)*2/steps;
maxSteps = steps*2;


% First approximation for first ray
[ip0,ip1,intersect] = circleIntersect(phasePos(1),s.P(1,:),s.V(1,:));
if intersect
    dv = (ip1-ip0);
    d = sqrt(dv(1)^2+dv(2)^2);
    n(2) = n0 + (phaseVal(2)/d);
else
    disp('Error: Ray did not intersect volume')
    return
end

% First approximation for the remaining rays
for rayInd = 2:s.nRays
    n(rayInd+1) = firstApproximation(rayInd);
end
nf = n;
s.resetSource();

% Second approximation for all rays exept last ray
% for rayInd = 1:s.nRays-1 % Exclude last ray that needs a special version
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
        dPhase(rayInd) = phaseVal(rayInd+1) - sum(dPhase);
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
            [p0,~,~] = circleIntersect(phasePos(1),s.P(rayInd,:),s.V(rayInd,:));
            rayPath(1,:) = p0;              % Set origin point
            v0 = s.V(rayInd,:);             % Current Vector
            rayStep = 2;                    % Current ray step
            minR = inf;                     % Current smalest radius
            cShell = 1;                     % Current shell nr.
%             r0 = sqrt(p0(1)^2 + p0(2)^2);   % Current distance from center.
            % Trace the ray through the volume until exit or max steps are
            % reached
%             shellPath = zeros(maxSteps,2).*nan;
%             shellPath(1,:) = p0;
            totalPath = 0;
            phaseSum_test = 0;
            
            % Initiate trace loop
            exitVolume = 0;     % Exit flag
            while ~exitVolume && rayStep < maxSteps
                % Determine Runge-Kutta matrix for current shell
                k = (n(rayInd) - n(rayInd+1)) / (shellR(rayInd) - shellR(rayInd+1));
                m = n(rayInd+1) - k*shellR(rayInd+1);
                nFunc = @(r) k*sqrt(r(1)^2 + r(2)^2) + m;
                D = @(r) k*(k + m/sqrt(r(1)^2 + r(2)^2)).*[r(1) r(2)];                                               
                
                % trace though current shell
                exitShell = 0;
                while ~exitShell && rayStep < maxSteps                
                    v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
                    A = ds.*D(p0);
                    B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
                    C = ds.*D(p0 + ds.*v0 + ds/2.*B);
                    p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
                    v1 = v0 + 1/6.*(A + 4.*B + C);
%                     if (norm(p1-p0) - ds)/ds > 10^-3
%                         disp('Error: steplength error')
%                         return
%                     end
                    if v1(2) > 0
                        phasePath = inf;
                        exitVolume = 1;
                        break
                    end
                    v1 = v1./sqrt(v1(1)^2 + v1(2)^2);
                    r1 = sqrt(p1(1)^2 + p1(2)^2);  % Current distance from center of volume
                    rayPath(rayStep,:) = p1;
                    rayStep = rayStep + 1;
                    dt = sqrt((p0(1)-p1(1))^2 + (p0(2)-p1(2))^2);
                    totalPath = totalPath + dt;
                    
                    if cShell == rayInd
                        phasePath = phasePath + dt;
                        phaseSum_test = phaseSum_test + dt*(nFunc(p0) - n(1));
%                         shellPath(rayStep,:) = p1;
                    else
                        phaseSum = phaseSum + dt*(nFunc(p0) - n(1));
                    end
                    v0 = v1;                         % Set new vector
                    p0 = p1;                         % Set new curent point
%                     r0 = r1;
                    if r1 < minR
                        minR = r1;   % New min distance to center of colume
                    end
                    if r1 <= shellR(cShell+1) && cShell < nShells
                        cShell = cShell + 1;
                        exitShell = 1;
                    elseif r1 > shellR(cShell)
                        cShell = cShell - 1;
                        exitShell = 1;
                        if cShell < 1
                            exitVolume = 1;
                        end
                    end
                end
            end
            figure(1)
            hold on; axis equal
            plotCircle(phasePos(1),2*pi,1)
            plot(rayPath(:,1),rayPath(:,2))
            
            % Calculate phase shift, back project ray & compare phase shift
            % values
            [phaseTot] = findTotPhase(p0,v0);
            ns  = n(1) + (phaseTot - phaseSum)/phasePath;
            n1 = (n(1+rayInd) + ns)/2;
            disp(strcat('Itteration: ',num2str(iter),', n1 = ', num2str(n1)))
            n(1 + rayInd) = n1;
            accuarcy = abs(nRef - n(1+rayInd));
            nRef = n(1+rayInd);
            shellR(1+rayInd) = minR;
            iter = iter + 1;
        end
        figure(2)
        hold on; axis equal
        plot(rayPath(:,1),rayPath(:,2))
        plotCircle(shellR(1),pi,2)
        plotCircle(shellR(rayInd+1),pi,2)
        plotLine([0 0],[shellR(1) 0],'k',2)
        plotLine([0 shellR(1)],[0 -shellR(1)],'k',2)
        
        s.PT{rayInd} = rayPath;
        s.totalPath(rayInd) = totalPath;
        s.phase(rayInd) = phaseSum_test;
    end

    function [phaseTot] = findTotPhase(p0,v0)
        t = -p0(2)/v0(2);
        ip = p0(1) + t*v0(1); % x-axis intersection distance        
        % find colsest phaseShift values and interpolate
        indP = find((phasePos-ip) > 0,1,'last');
        if isempty(indP)
            indP = 1;
        elseif indP == length(phasePos)
            indP = length(phasePos)-1;
        end
        k = (phaseVal(indP+1) - phaseVal(indP+1)) / (phasePos(indP+1) - phasePos(indP));
        phaseTot = phaseVal(indP) + k*(ip - phasePos(indP)); 
    end
end













