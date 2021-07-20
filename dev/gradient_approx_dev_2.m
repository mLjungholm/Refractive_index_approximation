clc
close all

lambda = 550*10^-9;
r = lambda*12*2;
n0 = 1.3; 
n1 = 1.45;
k = (n0-n1)/r^2;
nFunc = @(x) k.*x.^2 + n1; %<- True n(r) (parabolic)

% Linear refractive indices in shell 1
rayInd = 2;
ra = mPhasePos(rayInd);
rb = mPhasePos(rayInd+1);
na = nFunc(ra);
nb = relativeN(ra,rb,na,nFunc(rb),0);
lFunc = @(x) (na-nFunc(rb))/(ra-rb).*(x - ra) + na;

% sTrue = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);
% sLinearContr = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);
% sApproxTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);
sApproxTest2 = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);
steps = 10^3;

% ray_trace(sTrue,steps,n0,n1,r,'parabolic','meggit')
% ray_trace(sLinearContr,steps,n0,nb,r,'linear','meggit')

figure(1)
hold on; axis equal; grid on
sTrue.plotRays('--b',rayInd)
sLinearContr.plotRays('--r',rayInd)

% plotCircle(r,2*pi)
% plotCircle(mPhasePos(rayInd + 1),2*pi)

% plotLine([-r 0].*1.2,[r 0]*1.2,'k')
% plotLine([ 0 -r].*1.2,[0 r]*1.2,'k')

% sApproxTest.id = 'sApproxTest';
% nh = homogenious_approximation(sApproxTest, mPhaseShift, mPhasePos,n0); 
% shellR = [];
% [n,~,shellR,~] = gradient_approximation(sApproxTest,r,mPhasePos,mPhaseShift,nh,steps,sLinearContr,lFunc(rb),1,shellR);
[n2,~,shellR2,itterations] = gradient_approximation(sApproxTest2,r,mPhasePos,mPhaseShift,n,steps,sLinearContr,lFunc(shellR(3)),2,shellR);
% sApproxTest.plotRays('m',1)
% plotCircle(shellR(1 + 1),2*pi,[0.1,0.6,0.6])

% % Refractive index linearization
% figure(2)
% hold on; grid on
% x = linspace(ra,rb,100);
% plot(x,nFunc(x),'b')
% lFunc = @(x) (na-nFunc(rb))/(ra-rb).*(x - ra) + na;
% plot(x,lFunc(x),'r')

% % Phase shift profiles
% figure(3)
% hold on; grid on
% plot(truePhasePos,truePhaseShift,'b')
% plot(mPhasePos,mPhaseShift,'r')
% legend('True phase shift','Measured phase shift')

% Refractive index plot
% figure(4)
% hold on; grid on
% x = linspace(r,0,100);
% plot(x,nFunc(x),'b')
% scatter(shellR,n,'ro')

iterText = string(1:1:itterations(2));
legendText = ["Parabolic","Linear",iterText];
legend(legendText);

function [n,shellGradient,shellR,itterations] = gradient_approximation(source,r_edge,mPhasePos,mPhaseShift,n_homogenious,steps,controlSource,nContr,rayInd,shellR)
% n = zeros(source.nRays+1,1);
% n(1) = n0;
n = n_homogenious;
if isempty(shellR)
    shellR = mPhasePos;
end
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
        maxIter = 100;      % Maximum number of itterations
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
            phase_in_shellS = 0;
            % Initiate trace loop
            exitVolume = 0;     % Exit flag
            while ~exitVolume && rayStep < maxSteps
                % Find local curvature of radius
                v0 = v0./sqrt(v0(1)^2 + v0(2)^2);
                [p1,v1] = intStep(p0,v0,cShell);
                if isnan(p1(1))
                    stopflag = 1;
                end
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
                phase = phase + phaseShiftKnownGradient(r0,r1,cShell,ds);
                if cShell == rayInd
%                     phasePath = phasePath + dt;
                    phasePath = phasePath + ds;
                    shellPath(rayStep,:) = p1;
                    phase_in_shellS = phase_in_shellS + phaseShiftKnownGradient(r0,r1,cShell,ds);
                else
                    phaseSum = phaseSum + phaseShiftKnownGradient(r0,r1,cShell,ds);
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
            dPhase = phaseTot - phaseSum;
            if dPhase < 0
                n1 = n(rayInd+1);
                n2 = n1;
                nt = n1;
            else
                n1 = n(1) + (dPhase)/phasePath;
                nt = n1;
                n1 = getLinearGradient_CurvedTrace(shellR(rayInd),shellR(rayInd+1),n(rayInd),n1,dPhase,shellPath,ds);
                %             nR = relativeN(shellR(rayInd),shellR(rayInd+1),n(rayInd),n1,minR);
                n2 = (n(1+rayInd) + n1)/2;
            end
            
            if rayInd == 2
                
                fprintf('\n\n-------------------------------------------------\n')
                fprintf('Itteration : %u \n',iter);
                fprintf('Back Traced Phase  = %.4f \n',phaseTot*10^7);
                fprintf('Measured phase     = %.4f \n',phase*10^7);
                fprintf('Phase in shellS    = %.4f \n',phase_in_shellS*10^7);
                fprintf('Phase in outerS    = %.4f \n',phaseSum*10^7);
                fprintf('phaseT-phaseSum    = %.4f \n\n',dPhase*10^7);
                
                fprintf('Test source path in shell-S = %.4fum \n\n',totalPath*10^6);
                                                
                fprintf('n-Contr = %.4f \n',nContr);
                rN = relativeN(shellR(rayInd),minR,n(rayInd),n1,mPhasePos(rayInd+1));
                fprintf('n-Test  = %.4f \n',rN);
                nt = relativeN(shellR(rayInd),minR,n(rayInd),nt,mPhasePos(rayInd+1));
                fprintf('n-avrg  = %.4f \n',nt);
                fprintf('n-test/n-Contr = %.4f%% \n\n',(rN-nContr)/nContr*100)
                
                fprintf('R-min = %.4fum \n',minR*10^6)
                fprintf('R-s   = %.4fum \n\n',mPhasePos(rayInd+1)*10^6)
                
                plot(rayPath(:,1),rayPath(:,2))
            end
            
            n(1 + rayInd) = n2;
            accuarcy = abs(nRef - n(1+rayInd));
            nRef = n1;
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



