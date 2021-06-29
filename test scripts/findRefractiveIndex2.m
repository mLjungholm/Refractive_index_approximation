% Ray tracing based on algorithm by Nilsson, D-E et.al (1983), 
% using an approximation by Meggit & Meyer-Rochow (1975)
% Written by Mikael Ljungholm (2021)

function [n,shellR2] = findRefractiveIndex2(s,pp,ps,n0,steps)
n = zeros(s.nRays,1);
n(1) = n0;
shellR = pp(1:end-1);
shellR2 = shellR;
ds = pp(1)*2/steps;
maxSteps = steps*2;

% First approximation for first ray
% [ip0,ip1,intersect] = circleIntersect(pp(1),s.P(1,:),s.V(1,:));
% if intersect
%     dv = (ip1-ip0);
%     d = sqrt(dv(1)^2+dv(2)^2);
%     n(2) = n0 + (ps(2)/d);
% else
%     disp('Error: First aproximation intersect error')
%     return
% end

% First approximation for the remaining rays
% for rayInd = 2:s.nRays
%     n(rayInd+1) = firstApproximation(rayInd);
% end
% nf = n;
% s.resetSource();

% Second approximation for all rays exept last ray
% for rayInd = 1:s.nRays-1 % Exclude last ray that needs a special version
for rayInd = 1:4
    [ip0,~,intersect] = circleIntersect(pp(1),s.P(rayInd,:),s.V(rayInd,:));
    if intersect
        loopRay(ip0,rayInd);
    end
end

    % Calculates the first approximation of the refractive index ine each
    % shell. This assumes that the rays travel straight threw the slice.
%     function nr = firstApproximation(rayInd)
%         nShells = rayInd; % Number of shells to be taken into acount
%         dShells = zeros(nShells,1); % Distance traveled in each shell
%         for shellInd = 1:nShells
%             [ip0,ip1,~] = circleIntersect(shellR(shellInd),s.P(rayInd,:),s.V(rayInd,:));
%             dv = (ip1-ip0);
%             d = sqrt(dv(1)^2+dv(2)^2); % distance traveled in shell
%             dShells(shellInd) = d;
%         end
%          % Substract the distance from inner shells.
%         for shellInd = 1:nShells-1
%             dShells(shellInd) = dShells(shellInd) - dShells(shellInd+1); 
%         end
%         dPhase = zeros(nShells-1,1);
%         for shellInd = 1:nShells-1
%             dPhase(shellInd) = dShells(shellInd)*(n(1+shellInd) - n(1));
%         end
%         dPhase(rayInd) = ps(rayInd+1) - sum(dPhase);
%         nr = dPhase(end)/dShells(end) + n0;
%     end

    function loopRay(ip,rayInd)
%         close all
        maxIter = 10000;      % Maximum number of itterations
        threshold = 0.001;  % Threshold value
        iter = 1;           % Current itteration
        accuarcy = inf;     % Current accuarcy
        nRef = 0;           % Current calculated refractive index
        nShells = rayInd;   % Number of shells the ray needs to pass
        rayPath = [];       % Path of ray
        % Run untill threshold is met or max itterations are reached
        n(rayInd +1) = n(rayInd)*1.001;
        disp(strcat('rayNr : ',num2str(rayInd)))
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
%             k = (n(shellInd+1)-n(shellInd))/(shellR(shellInd +1) - shellR(shellInd));
%             phaseFunc = @(x) ds*(k*(x - shellR(shellInd)) + n(shellInd) - n(1));
            while ~exitVolume && rayStep < maxSteps
                % Find local curvature of radius
                v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
                
                dr = findLocalRadius(n(cShell),n(cShell+1),(shellR(cShell)-shellR(cShell+1)),p0,v0);
                if dr < 0
                    errflag = 1;
                end
                [p1,v1] = curvePath(p0,v0,dr,ds); % Move one integration distance
                if v1(2) > 0
                    errflag = 1;
                    minR = rayPath(1,1);
                    p0 = [shellR(1) 0];
                    v0 = [1,-1];
                    continue
                elseif v1(1) > v0(1)
                    errflag = 1;
                end
                
                r1 = sqrt(p1(1)^2 + p1(2)^2);  % Current distance from center of volume
                rayStep = rayStep + 1;
                rayPath(rayStep,:) = p1;         % Add step to ray path
                
                %                 if cShell == rayInd
                %                     phasePath = phasePath + ds;
                %                 else
                                    phaseSum = phaseSum + phaseShiftKnownGradient((r1+r0)/2,cShell);
                %                 end
%                 phaseSum = phaseSum + phase((r1+r0)/2,cShell);
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
            
             
            
            [phaseTot,xip] = findTotPhase(p0,v0);
            if phaseSum < phaseTot
                n1 = n(rayInd +1)*1.0005;
            else
                n1 = n(rayInd +1)*0.9995;
            end
%             ns  = n(1) + (phaseTot - phaseSum)/phasePath;           
%             n1 = getLinearGradient_CurvedTrace(shellR(rayInd),shellR(rayInd+1),n(rayInd),ns,(phaseTot - phaseSum),[rayPath(2:end,1) rayPath(2:end,2)],ds);
%             n1 = (n(1+rayInd) + n1)/2;
            
%             disp(strcat('Itteration: ',num2str(iter),', n1 = ', num2str(n1),' ns = ', num2str(ns), ' nOld = ',num2str(n(rayInd+1)),' r = ',num2str(minR))) 
%             disp('')
            n(1 + rayInd) = n1;
            accuarcy = abs(phaseTot - phaseSum)/phaseSum;
%             nRef = n(1+rayInd);
            shellR2(1+rayInd) = minR;
            iter = iter + 1;
            if rayInd == 4
                testflag = 1;
            end
        end
        figure(2)
        hold on; axis equal
        plot(rayPath(:,1),rayPath(:,2))
        plotCircle(shellR2(1),pi,2)
        plotCircle(shellR2(rayInd+1),pi,2)
        plotLine([0 0],[shellR2(1) 0],'k',2)
        plotLine([0 shellR2(1)],[0 -shellR2(1)],'k',2)
        
        s.PT{rayInd} = rayPath;
        disp(strcat('Itterations :',num2str(iter)))

        %         disp(strcat('RayNr: ',num2str(rayInd),', Itterations :',num2str(iter),', Accuarcy: ',num2str(accuarcy),', n = ',num2str(n(1 + rayInd))))
    end


    function r = findLocalRadius(n0,n1,d,p,v)
        sN = -p';
        sN = sN./sqrt(sN(1)^2 + sN(2)^2);
        theta = acos(sN(1)*v(1) + sN(2)*v(2));
        r = ((n0 + n1)/2) / (sin(theta)*((n1-n0)/d));
    end

    function phaseS = phaseShiftKnownGradient(r,shellInd)
        k = (n(shellInd+1)-n(shellInd))/(shellR(shellInd +1) - shellR(shellInd));        
        dn = k*(r - shellR(shellInd)) + n(shellInd);
        phaseS = ds*(dn-n(1));
    end

    function [k,m] = getSlope(n0,n1,r0,r1)
        k = (n1-n0)/(r1-r0);
        m = n0-(k*r0);
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













