% Development script for the gradient aproximation algorithm

% Load data or generate new interference data with the
% "createPhaseShift.m" script.

% create a source with variable spacing to match the positions of the
% measured phase values.
clc
close all

% Getting the max n1 value if the gradient was linear within the range of
% the first shell
k = (n0-n1)/r^2;
nFunc = @(x) k.*x.^2 + n1; %<- True n(r) (parabolic)

shellProfile = zeros(nPeaks - 1,3);
na = nFunc(r); %<- Shloud be na=n0 for fisrt shell
nb = nFunc(peakPos(1 + 1)); %<- max value in shell 1
ra = r;  
rb = peakPos(1 + 1);
kL = (na-nb)/(ra-rb); % Linear gradient of n(r) in shell 1 (dif(n(x))
nFuncL = @(x) kL.*(x-ra) + na; % Linear function in shell 1
nb = nFuncL(0); % Max value for lens if all the volume should have the same linear gradient as shell 1
shellProfile(1,:) = [na,nb,r];
for shellInd = 2:nPeaks - 1
na = nFunc(peakPos(shellInd)); %<- Shloud be na=n0 for fisrt shell
nb = nFunc(peakPos(shellInd + 1)); %<- max value in shell 1
ra = peakPos(shellInd);  
rb = peakPos(shellInd + 1);
kL = (na-nb)/(ra-rb); % Linear gradient of n(r) in shell 1 (dif(n(x))
nFuncL = @(x) kL.*(x-ra) + na; % Linear function in shell 1
nb = nFuncL(0); % Max value for lens if all the volume should have the same linear gradient as shell 1
shellProfile(shellInd,:) = [na,nb,ra];
end


sLinearContr = Source_2d([peakPos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);
sLinearContr.id = 'sLinearContr';
sApproxTest = Source_2d([peakPos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);
sApproxTest.id = 'sApproxTest';
nh = homogenious_approximation(sApproxTest, mPhaseShift, peakPos,n0);
steps = 10^4;

% Ray 1
% Trace the test ray using the linearization of n(x)-true in shell 1
ray_trace_single_ray(sLinearContr,steps,n0,n1,r,shellProfile,'meggitStepGradient',1); 
% Approximating the refractive index

% shellProfile = shellProfile(2,1);

[n,~,shellR,~] = gradient_approximation(sApproxTest,r,peakPos,mPhaseShift,nh,steps,sLinearContr,shellProfile(2,1),1);


% % Ray 2
% Trace the test ray using the linearization of n(x)-true in shell 1
% ray_trace_single_ray(sLinearContr,steps,na,nb,r,shellProfile,'meggitStepGradient',2); 
% % Approximating the refractive index
% [n,~,shellR, ~] = gradient_approximation(sApproxTest,shellR,mPhaseShift,n,steps,sLinearContr,shellProfile(3,1),2);



figure(1)
hold on; axis equal; grid on
for rayInd = 1:1
sContr.plotRays('--r',rayInd)
sApproxTest.plotRays('m',rayInd)
sLinearContr.plotRays('--b',rayInd)
plotCircle(shellR(rayInd+1),2*pi,1)
end
plotCircle(r,2*pi,1)
plotLine([-r 0].*1.2,[r 0]*1.2,'k')
plotLine([ 0 -r].*1.2,[0 r]*1.2,'k')
% iterText = string(1:1:itterations(1));
% legendText = [iterText "Control" "Linear",];
% legend(legendText);



% figure(2)
% hold on; grid on
% scatter(shellR,n,'bo')
% plot(linspace(r,0,100),nFunc(linspace(r,0,100)),'k')
% for i = 1:size(shellProfile,1)-1
%     x = [shellProfile(i,3);shellProfile(i+1,3)];
%     y = [shellProfile(i,1);shellProfile(i+1,1)];
%     plot(x,y,'r');
% end
% x = [shellProfile(end,3);0];
% y = [shellProfile(end,1);shellProfile(end,2)];
% plot(x,y,'r');

% Ray tracing based on algorithm by Nilsson, D-E et.al (1983),
% using an approximation by Meggit & Meyer-Rochow (1975)
% Written by Mikael Ljungholm (2021)

function [n,shellGradient,shellR,itterations] = gradient_approximation(source,r_edge,mPhasePos,mPhaseShift,n_homogenious,steps,controlSource,nContr,rayInd)
% n = zeros(source.nRays+1,1);
% n(1) = n0;
n = n_homogenious;
shellR = mPhasePos;
ds = mPhasePos(1)*2/steps;
maxSteps = steps*2;
shellGradient = zeros(source.nRays,2).*nan;
itterations = zeros(source.nRays,1);

setGradient(1);
% Second approximation for all rays exept last ray
% for rayInd = 1:s.nRays-1 % Exclude last ray that needs a special version
% for rayInd = 1:s.nRays-1
% for rayInd = 1:source.nRays
% for rayInd = 1:2
   [~] = loopRay(rayInd);
%    if stopflag
%        break
%    end
% end


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
                    if cShell > 1
                    cShell = cShell - 1;
                    elseif r1 > r_edge
%                     if cShell < 1
                        exitVolume = 1;
                    end
                    
                end
            end
            
            % Calculate phase shift, back project ray & compare phase shift
            % values
            phaseTot = findTotPhase(p0,v0,mPhasePos,mPhaseShift);
            n1  = n(1) + (phaseTot - phaseSum)/phasePath;
%             nR = getLinearGradient_CurvedTrace(shellR(rayInd),minR,n(rayInd),n1,phaseTot,rayPath,ds);
            nR = relativeN(shellR(rayInd),shellR(rayInd+1),n(rayInd),n1,minR);
            n2 = (n(1+rayInd) + nR)/2;
            if rayInd == 1

            % Code for pringing values in each itteration
            fprintf('\n\nItteration : %u \n',iter);
            
%             fprintf('minR/shellR = %.4f%% \n',(minR-mPhasePos(rayInd+1))/mPhasePos(rayInd+1)*100)
            
            fprintf('Back Traced Phase T-testS: %.4f \n',phaseTot*10^7);
            T_contr = findTotPhase(controlSource.P(rayInd,:),controlSource.V(rayInd,:),mPhasePos,mPhaseShift);
            fprintf('Back Traced Phase T-contr: %.4f \n',T_contr*10^7);
            fprintf('Measured Contr Phase T-m:  %.4f \n\n',controlSource.phase(rayInd)*10^7);
            
            fprintf('Test source path Length   d-testS: %.4f \n',totalPath*10^7);
            fprintf('Contr source path Length  d-testS: %.4f \n\n',controlSource.totalPath(rayInd)*10^7);
            
            fprintf('n-True: %.4f \n',nContr);
            rN = relativeN(shellR(rayInd),minR,n(rayInd),nR,mPhasePos(rayInd+1));
            fprintf('n-TestM: %.4f \n',rN);
            fprintf('n-TestM: %.4f for r = %.2fum\n',nR,minR*10^6);
            na = T_contr/controlSource.totalPath(rayInd) + n(1);
            [~,d_contr] = getClosestPoint(controlSource,rayInd,[0 0]);
            rN = relativeN(shellR(rayInd),d_contr,n(rayInd),na,mPhasePos(rayInd+1));
            fprintf('n-contr: %.4f for r = %.2fum -> n = %.4f \n',na,d_contr*10^6,rN);
            fprintf('n-test/n-True= %.4f%% \n',(rN-nContr)/nContr*100) 
            end
            
            n(1 + rayInd) = n2;
            accuarcy = abs(nRef - n(1+rayInd));
            nRef = nR;
            shellR(1+rayInd) = minR;
            setGradient(rayInd);
            iter = iter + 1;
            if iter == maxIter
                stopflag = 1;
                fprintf('Error in gradient_approximation()\n')
                fprintf('Exceeding maximum itterations! \n Faulty ray ind = %u \n\n',rayInd)
                break
            end            
        end
        if n(1+rayInd) < 0
            fprintf('Error: approximated refractive index was negative! \n Faulty rayInd = %u \n',rayInd)
            n(1+rayInd) = n_homogenious(1+rayInd);
            stopflag = 1;
            return
        end
%         shellR(1+rayInd) = minR;
%         nR = relativeN(shellR(rayInd),mPhasePos(rayInd+1),n(rayInd),n(rayInd+1),minR);
%         nR = relativeN(shellR(rayInd),minR,n(rayInd),n(rayInd+1),mPhasePos(rayInd+1));
%         n(rayInd+1) = nR;
%         setGradient(rayInd);
        itterations(rayInd) = iter-1;
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

%     function phaseTot = findTotPhase(p0,v0)
%         t = -p0(2)/v0(2);
%         xip = p0(1) + t*v0(1); % x-axis intersection distance
%         % find colsest phaseShift values and interpolate
%         pind = find((mPhasePos-xip)>0,1,'last');
%         if isempty(pind)
%             pind = 1;
%         elseif pind == length(mPhasePos)
%             pind = length(mPhasePos) - 1;
%         end
%         p0 = mPhasePos(pind);
%         p1 = mPhasePos(pind+1);
%         shift0 = mPhaseShift(pind);
%         shift1 = mPhaseShift(pind+1);       
%         phaseTot = (shift1-shift0)/(p1-p0)*(xip-p0) + shift0;
%         if phaseTot < 0
%             phaseTot = 0;
%         end
%     end
end

    function phaseTot = findTotPhase(p0,v0,mPhasePos,mPhaseShift)
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

    
    function rN = relativeN(r0,r1,n0,n1,r)
    kL = (n0-n1)/(r0-r1); % Linear gradient of n(r) in shell 1 (dif(n(x))
    nFuncL = @(x) kL.*(x-r0) + n0; % Linear function in shell 1
    rN = nFuncL(r);
    end


    function [closestPoint,d] = getClosestPoint(source,rayNr,point)
    path = source.path{rayNr};
    key = ~isnan(path(:,1));
    path = path(key,:);
    steps = size(path,1);
    d = inf;
    closestPoint = [nan nan];
    for ind = 1:steps
        testD = norm(path(ind,:)-point);
        if testD < d
            d = testD;
            closestPoint = path(ind,:);
        end
    end
    end




