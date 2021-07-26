clc
close all

%-------------------------------------------------------------------------%
                            % True gradient %
lambda = 550*10^-9;
r = lambda*12*2;
n0 = 1.3; 
n1 = 1.45;
k = (n0-n1)/r^2;
nFunc = @(x) k.*x.^2 + n1; %<- True n(r) (parabolic)

%-------------------------------------------------------------------------%

% Test source
sTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% Control source
% sLinear = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% Slice object
C = IndexSlice(mPhasePos,mPhaseShift,n0);

% Homogenious approximation
homogenious_approximation(sTest,C); 

rayInd = 1;
steps = 10^3;
analytical = gradient_approximation(sTest, C, steps, rayInd, nFunc, sContr, 0);
% figure(4)
% plot(na)
rayInd = 2;
analytical = gradient_approximation(sTest, C, steps, rayInd, nFunc, sContr,0);

% figure(1)
% hold on; axis equal; grid on
% sTest.plotRays('r',1)
% sContr.plotRays('--b',1)
% sTest.plotRays('r',2)
% sContr.plotRays('--b',2)
% sLinear.plotRays('--g',1)

% figure(2)
% hold on; grid on
% C.plotNCurve
% C.plotN(1)
% plot(linspace(r,0,100),nFunc(linspace(r,0,100)),'b')

figure(3)
hold on; grid on

x = 1:length(analytical{1});
plot(x,analytical{1},x,analytical{2},x,analytical{3},x,analytical{4},x,analytical{5})
legend("ns","n1","n2","nTrue","n0")

% iterText = string(1:1:itterations(2));
% legendText = ["Parabolic","Linear",iterText];
% legend(legendText);

function analytical = gradient_approximation(source, C, steps, rayInd, nTrue, sContr, print_Analytical)

ds = C.r(1)*2/steps;
maxSteps = steps*2;
itterations = zeros(source.nRays,1);

C.calculateGradient('all');

[~] = loopRay(rayInd);


    function stopflag = loopRay(rayInd)
        stopflag = 0;
        % Set up recurcive aproximation
        maxIter = 30;      % Maximum number of itterations
        threshold = 0.001;  % Threshold value
        
        iter = 1;           % Current itteration
        accuarcy = inf;     % Current accuarcy
        nRef = C.n(rayInd+1); % Current calculated refractive index
        
        nShells = rayInd;   % Number of shells the ray needs to pass
        rayPath = [];       % Path of ray
        
        analytical = cell(4,1); % Test Var
        Ans = zeros(maxIter,1).*nan; % ns vals
        An1 = zeros(maxIter,1).*nan; % n1 vals
        An2 = zeros(maxIter,1).*nan; % n2 vals
        AnT = zeros(maxIter,1).*nan; % n-True vals        
        Anb = zeros(maxIter,1).*nan; % n(-1) vals
        Atp = zeros(maxIter,1).*nan; % totalPhase vals
        Abp = zeros(maxIter,1).*nan; % back phase vals
        Aip = zeros(maxIter,1).*nan; % inner phase vals
        Aop = zeros(maxIter,1).*nan; % outer phase vals
        Aac = zeros(maxIter,1).*nan; % accuarcy vals
        
        % Run untill threshold is met or max itterations are reached        
        while iter < maxIter && accuarcy > threshold
            
            rayPath = zeros(maxSteps,2).*nan;  % Create new path vector            
            testShell_d = 0;        % Total distance in shell to be calculated
            outerPhaseSum = 0;      % Sum of phase shift in outer shells
            minR = inf;               % Smalest distance to center for the path
            rayStep = 1;            % Current ray step
            cShell = 1;             % Current shell nr.
            
            % Extra test variables
            testShellPath = zeros(maxSteps,2).*nan; % Path in the test shell
            totalPath = 0;       % Total path in all shells
            totalPhase = 0;      % Total phase shift in all shells
            testShellPhase = 0;  % Phase shift in test shell
            outerShell_d = 0;    % Total path in outer shells
             
            
            % Intersection with outer shell
            [p0,~,~] = circleIntersect(C.r(1),source.P(rayInd,:),source.V(rayInd,:));
            
            rayPath(rayStep,:) = p0;              % Set origin point
            v0 = source.V(rayInd,:);        % Current Vector
            r0 = norm(p0);   % Current distance from center.
                                                            
            % Trace the ray through the volume until exit or max steps are
            % reached            
            % Initiate trace loop
            exitVolume = 0;     % Exit flag
            while ~exitVolume && rayStep < maxSteps
                rayStep = rayStep + 1;
                
                % Find local curvature of radius
                v0 = v0./sqrt(v0(1)^2 + v0(2)^2); % Doulbe cheking normalization of direction vector.
                [p1,v1] = intStep(p0,v0,cShell);  % One integration step
                
                if v1(2) > 0 % If v1(2) > 0 -> ray has looped around (to high n-index)
                    minR = rayPath(1,1);
                    p0 = [shellR(1) 0];
                    v0 = [1,-1];
                    continue
                end
                
                r1 = sqrt(p1(1)^2 + p1(2)^2);     % Current distance from center of volume              
                rayPath(rayStep,:) = p1;          % Add step to ray path
                totalPath = totalPath + ds;       % Add step-length to total path (test var)
                dP = phaseShift(r0,r1,cShell);    % Calculate phase shift for step
                totalPhase = totalPhase + dP;     % Add the phase shift for the step to the total phase (test var)
                
                if cShell == rayInd                         % If the current shell is the test shell
                    testShell_d = testShell_d + ds;         % Add step length to path in test shell  
                    testShellPath(rayStep,:) = p1;          % Add step point path in shell (test var)        
                    testShellPhase = testShellPhase + dP;   % Add phase shift to total phase in test shell (test var)                    
                else
                    outerPhaseSum = outerPhaseSum + dP;     % Add phase shift to outer shell phase sum
                    outerShell_d = outerShell_d + ds;       % Add step length to outer shell path (test var)
                end
                
                v0 = v1;    % Set new current vector
                p0 = p1;    % Set new curent point
                r0 = r1;    % Set new current r
                
                if r1 < minR % Check if new min distance to center of volume
                    minR = r1;
                end
                
                
                % Test if to move to lower shell or exit volume                
                if r1 <= C.r(cShell+1) && cShell < nShells                    
                    cShell = cShell + 1; % Go down one shell
                
                elseif r1 > C.r(cShell)  % r lager than current shell                    
                    if cShell > 1        % Not in outer most shell
                    cShell = cShell - 1; % Go up one shell
                    elseif r1 > C.r(1)   % r larger than outer most shell
                        exitVolume = 1;  % Exit volume
                    end
                end
            end  % End of trace loop -------------------------------------%
            
            %-------------------------------------------------------------%
            % Calculate phase shift, back project ray & compare phase shift
            % values.
            
            rIp = projectRay(p0,v0);                 % Find ray back-projection to focus plane           
            projectedPhase = C.getPhaseShift(rIp);   % Get phase for the backprojected radius            
            
            dPhase = (projectedPhase - outerPhaseSum);
            if dPhase < 0
                dPhase = 0;            
            end
            if dPhase == 0
                ns = (nRef + C.n(1))/2;
            else
                ns = dPhase/testShell_d + C.n(1);
            end
            n1 = getLinearGradient_CurvedTrace(C.r(rayInd),C.r(rayInd+1),C.n(rayInd),ns,dPhase,testShellPath,ds);              
%             n1 = getLinearGradient_CurvedTrace(C.r(rayInd),minR,C.n(rayInd),ns,dPhase,testShellPath,ds); 
            n2 = (C.n(rayInd +1) + n1)/2; % Create an average n-index of the last two
            if n2 < 1
                fprintf('Error: approximated refractive index was less than 1! \n Faulty rayInd = %u \n',rayInd)
                break
            end
            C.n(1 + rayInd) = n2;
            C.calculateGradient(rayInd);
            
            accuarcy = abs(nRef - n1);
            nRef = n1;
            if print_Analytical
                printVals(iter, projectedPhase, totalPhase, testShellPhase, outerPhaseSum, dPhase, totalPath, rayPath, minR, n1)
            end
            Ans(iter) = ns;
            An1(iter) = n1;
            An2(iter) = n2;
            AnT(iter) = nTrue(C.r(rayInd+1));
            Anb(iter) = C.n(rayInd);
            Atp(iter) = totalPhase;
            Abp(iter) = projectedPhase;
            Aip(iter) = testShellPhase;
            Aop(iter) = outerPhaseSum;
            Aac(iter) = accuarcy;                        
            iter = iter + 1;
            if iter == maxIter
                stopflag = 1;
                fprintf('Error in gradient_approximation()\n')
                fprintf('Exceeding maximum itterations! \n Faulty ray ind = %u \n\n',rayInd)
                break
            end            
        end % End of Itteration loop -------------------------------------%
        
        C.n(rayInd+1) = C.getN(rayInd,minR);
        C.r(rayInd+1) = minR;
        C.calculateGradient(rayInd);
        C.calculateGradient(rayInd+1);
        
        analytical = {Ans, An1, An2, AnT, Anb, Atp, Abp, Aip, Aop, Aac, iter-1};
        
        itterations(rayInd) = iter-1;
        source.P(rayInd,:) = p1;
        source.V(rayInd,:) = v1;
        source.path{rayInd} = rayPath;
        source.phase(rayInd) = totalPhase;
        source.totalPath(rayInd) = totalPath;
    end

    function [p1,v1] = intStep(p0,v0,shellInd)
        r0 = norm(p0);
        gV = -p0./r0;
        gV = gV./norm(gV);
        theta = acos(dot(v0,gV)./norm(v0));
        nr = C.getN(shellInd,r0);
        r = nr/(sin(theta) * -C.gradient(shellInd));
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

    function dP = phaseShift(r0,r1,shellInd)
        nr = C.getN(shellInd,(r0+r1)/2);
        dP = ds*(nr-C.n(1));
    end

    function rIp = projectRay(p,v)
        t = -p(2)/v(2);
        rIp = p(1) + t*v(1);
    end

    function printVals(iter, projectedPhase, totalPhase, testShellPhase, outerPhaseSum, dPhase, totalPath, rayPath, minR, n1)
        fprintf('\n\n-------------------------------------------------\n')
        fprintf('Itteration : %u \n',iter);
        fprintf('Back Traced Phase  = %.4f \n',projectedPhase*10^7);
        fprintf('Measured phase     = %.4f \n',totalPhase*10^7);
        fprintf('Phase in shellS    = %.4f \n',testShellPhase*10^7);
        fprintf('Phase in outerS    = %.4f \n',outerPhaseSum*10^7);
        fprintf('phaseT-phaseSum    = %.4f \n',dPhase*10^7);
        fprintf('(parabolic) trace  = %.4f \n\n',sContr.phase(rayInd)*10^7);
        
        fprintf('Test source path in shell-S = %.4fum \n\n',totalPath*10^6);
        
        fprintf('n-Contr = %.4f \n',nTrue(C.r(rayInd+1)));
%         rN = relativeN(C.r(rayInd),minR,C.n(rayInd),n1,C.r(rayInd+1));
        fprintf('n-Test  = %.4f \n',n1);
        fprintf('n-test/n-Contr = %.4f%% \n\n',(n1-nTrue(C.r(rayInd)))/nTrue(C.r(rayInd))*100)
        
        fprintf('R-min = %.4fum \n',minR*10^6)
        fprintf('R-s   = %.4fum \n\n',C.r(rayInd+1)*10^6)
        
%         figure(9)
%         hold on; axis equal; grid on
%         plot(rayPath(:,1),rayPath(:,2))        
    end

end

%-------------------------------------------------------------------------%
                    % External functions %

                    
                   
                    
                    
% function phaseTot = findTotPhase(p0,v0,mPhasePos,mPhaseShift)
% t = -p0(2)/v0(2);
% xip = p0(1) + t*v0(1); % x-axis intersection distance
% % find colsest phaseShift values and interpolate
% pind = find((mPhasePos-xip)>0,1,'last');
% if isempty(pind)
%     pind = 1;
% elseif pind == length(mPhasePos)
%     pind = length(mPhasePos) - 1;
% end
% p0 = mPhasePos(pind);
% p1 = mPhasePos(pind+1);
% shift0 = mPhaseShift(pind);
% shift1 = mPhaseShift(pind+1);
% phaseTot = (shift1-shift0)/(p1-p0)*(xip-p0) + shift0;
% if phaseTot < 0
%     phaseTot = 0;
% end
% end
% 
% 
function rN = relativeN(r0,r1,n0,n1,r)
kL = (n0-n1)/(r0-r1); % Linear gradient of n(r) in shell 1 (dif(n(x))
nFuncL = @(x) kL.*(x-r0) + n0; % Linear function in shell 1
rN = nFuncL(r);
end
% 
% 
% function [closestPoint,d] = getClosestPoint(source,rayNr,point)
% path = source.path{rayNr};
% key = ~isnan(path(:,1));
% path = path(key,:);
% steps = size(path,1);
% d = inf;
% closestPoint = [nan nan];
% for ind = 1:steps
%     testD = norm(path(ind,:)-point);
%     if testD < d
%         d = testD;
%         closestPoint = path(ind,:);
%     end
% end
% end

 
% fprintf('\n\n-------------------------------------------------\n')
% fprintf('Itteration : %u \n',iter);
% fprintf('Back Traced Phase  = %.4f \n',projectedPhase*10^7);
% fprintf('Measured phase     = %.4f \n',totalPhase*10^7);
% fprintf('Phase in shellS    = %.4f \n',testShellPhase*10^7);
% fprintf('Phase in outerS    = %.4f \n',outerPhaseSum*10^7);
% fprintf('phaseT-phaseSum    = %.4f \n\n',dPhase*10^7);
% 
% fprintf('Test source path in shell-S = %.4fum \n\n',totalPath*10^6);
% 
% fprintf('n-Contr = %.4f \n',nContr);
% rN = relativeN(shellR(rayInd),minR,n(rayInd),n1,mPhasePos(rayInd+1));
% fprintf('n-Test  = %.4f \n',rN);
% nt = relativeN(shellR(rayInd),minR,n(rayInd),nt,mPhasePos(rayInd+1));
% fprintf('n-avrg  = %.4f \n',nt);
% fprintf('n-test/n-Contr = %.4f%% \n\n',(rN-nContr)/nContr*100)
% 
% fprintf('R-min = %.4fum \n',minR*10^6)
% fprintf('R-s   = %.4fum \n\n',mPhasePos(rayInd+1)*10^6)
% 
% plot(rayPath(:,1),rayPath(:,2))


