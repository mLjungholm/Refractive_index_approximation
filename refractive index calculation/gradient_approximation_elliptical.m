% In use
function gradient_approximation_elliptical(source, C, steps, varargin)

ds = C.r(1)*2/steps; % Set step size
maxSteps = steps*2;  % Set max steps

C.calculateGradient('all'); % Set the initial gradient for all shells

% Check if all rays or a single defined ray is to be traced.
if isequal(varargin{1},'all') || isempty(varargin)
    rayRange = 1:source.nRays;
else
    rayRange = varargin{1}:varargin{1};
end

% Main loop. sadly this canot be parallizised since the result of each ray
% builds on the previos ones.
for rayInd = rayRange
[~] = loopRay(rayInd);
end

    function stopflag = loopRay(rayInd)
        stopflag = 0;
        % Set up recurcive aproximation
        maxIter = 50;      % Maximum number of itterations
        threshold = 0.001;  % Threshold value        
        iter = 1;           % Current itteration
        
        nShells = rayInd;   % Number of shells the ray needs to pass
        rayPath = [];       % Path of ray
        
        Ans = zeros(maxIter,1).*nan; % ns vals
        AnH = zeros(maxIter,1).*nan; % n_High
        AnL = zeros(maxIter,1).*nan; % n_Low
        Atp = zeros(maxIter,1).*nan; % totalPhase vals
        Abp = zeros(maxIter,1).*nan; % back phase vals
        Aip = zeros(maxIter,1).*nan; % inner phase vals
        Aop = zeros(maxIter,1).*nan; % outer phase vals
        Aac = zeros(maxIter,1).*nan; % accuarcy vals
        
        % Run untill threshold is met or max itterations are reached        
        while iter < maxIter
            
            rayPath = zeros(maxSteps,2).*nan;  % Create new path vector            
            inner_d = 0;      % Total distance in shell to be calculated
            outerPhase = 0;   % Sum of phase shift in outer shells
            minR = inf;       % Smalest distance to center for the path
            rayStep = 1;      % Current ray step
            cShell = 1;       % Current shell nr.
            
            % Extra test variables
            testShellPath = zeros(maxSteps,2).*nan; % Path in the test shell
            totalPath = 0;    % Total path in all shells
            totalPhase = 0;   % Total phase shift in all shells
            innerPhase = 0;   % Phase shift in test shell
            outer_d = 0;      % Total path in outer shells
             
            
            % Intersection with outer shell
            [p0,~,~] = ellipseIntersect(C.r(1),C.r(1)*C.rMinor,source.P(rayInd,:),source.V(rayInd,:));
            
            rayPath(rayStep,:) = p0;    % Set origin point
            v0 = source.V(rayInd,:);    % Current Vector
            r0 = norm(p0);              % Current distance from center.
                                                            
            % Trace the ray through the volume until exit or max steps are
            % reached            
            % Initiate trace loop
            exitVolume = 0;     % Exit flag
            while ~exitVolume && rayStep < maxSteps
                rayStep = rayStep + 1;
                
                v0 = v0./sqrt(v0(1)^2 + v0(2)^2);      % Doulbe cheking normalization of direction vector.
                if iter == 1 || rayInd == source.nRays % Do not curve if it is the first itteration
                    p1 = p0 + ds.*v0;
                    v1 = v0;
                else
                    [p1,v1] = intStep(p0,v0,cShell);  % One curved integration step
                end
                
                if v1(2) > 0 % If v1(2) > 0 -> ray has looped around (to high n-index)
                    minR = rayPath(1,1);
                    p0 = [C.r(rayInd) 0];
                    v0 = [1,-1];
                    continue
                end
                
                r1 = sqrt(p1(1)^2 + p1(2)^2);     % Current distance from center of volume              
                rayPath(rayStep,:) = p1;          % Add step to ray path
                totalPath = totalPath + ds;       % Add step-length to total path (test var)
                dP = phaseShift(p0,cShell);    % Calculate phase shift for step
                totalPhase = totalPhase + dP;     % Add the phase shift for the step to the total phase (test var)
                
                if cShell == rayInd                 % If the current shell is the test shell
                    inner_d = inner_d + ds;         % Add step length to path in test shell  
                    testShellPath(rayStep,:) = p1;  % Add step point path in shell (test var)        
                    innerPhase = innerPhase + dP;   % Add phase shift to total phase in test shell (test var)                    
                else
                    outerPhase = outerPhase + dP;   % Add phase shift to outer shell phase sum
                    outer_d = outer_d + ds;         % Add step length to outer shell path (test var)
                end
                
                v0 = v1;    % Set new current vector
                p0 = p1;    % Set new curent point
                r0 = r1;    % Set new current r
                
                relativeR = C.getRelativeRadius(p1);
                if relativeR < minR % Check if new min distance to center of volume
                    minR = relativeR;
                end
                                
                % Test if to move to lower shell or exit volume                
                if relativeR <= C.r(cShell+1) && cShell < nShells                    
                    cShell = cShell + 1; % Go down one shell
                
                elseif relativeR > C.r(cShell)  % r lager than current shell                    
                    if cShell > 1        % Not in outer most shell
                    cShell = cShell - 1; % Go up one shell
                    elseif relativeR > C.r(1)   % r larger than outer most shell
                        exitVolume = 1;  % Exit volume
                    end
                end
            end  % End of trace loop -------------------------------------%
            
            %-------------------------------------------------------------%
            % Calculate phase shift, back project ray & compare phase shift
            % values.
            
            rIp = projectRay(p0,v0);                 % Find ray back-projection to focus plane           
            projectedPhase = C.getPhaseShift(rIp);   % Get phase for the backprojected radius            
            
            % Determine difference in measured phase and projected phase
            dPhase = (totalPhase - projectedPhase);
            
            % If it is the first trace determine max- and min-n.
            if iter == 1
                nLow = C.n(rayInd);
                nHigh = (projectedPhase*1.5 - outerPhase)/inner_d + C.n(1);
%                 nHigh = C.n(rayInd + 2);
                iter = iter + 1;
                continue
            end
            
            % Check accuarcy
            accuarcy = abs(dPhase/projectedPhase);
            if accuarcy < threshold
                break
            end
            
            % if dPhase > 0 -> reduce ns
            % if dPhase < 0 -> increase ns
            if dPhase > 0
                ns = (C.n(rayInd+1) + nLow)/2;
                nHigh = C.n(rayInd+1);                
            else
                ns = (C.n(rayInd+1) + nHigh)/2;
                nLow = C.n(rayInd+1);
            end
            
            C.n(rayInd + 1) = ns;        % Set new approximation
            C.calculateGradient(rayInd); % Uppdate gradient
                                    
            Ans(iter) = ns;
            AnH(iter) = nHigh;
            AnL(iter) = nLow;
            Atp(iter) = totalPhase;
            Abp(iter) = projectedPhase;
            Aip(iter) = innerPhase;
            Aop(iter) = outerPhase;
            Aac(iter) = accuarcy;  
            iter = iter + 1;   % Uppdate itteration
            if iter == maxIter % Check if exceding maximum itterations
                stopflag = 1;
                fprintf('Error in gradient_approximation()\n\n')
                fprintf('Ray exceeding maximum itterations! \nFaulty ray ind = %u \n',rayInd)
                fprintf('Accuarcy ended at %.2f%% \n\n',(1-accuarcy)*100)
                break
            end            
        end % End of Itteration loop -------------------------------------%
        
        % Set new radius and corresponding refractive index to the shell
        % (except for last ray. i.e ray passing through center)
        if rayInd < source.nRays
        C.n(rayInd+1) = C.getN(rayInd,minR); % Uppdate n-index
        C.r(rayInd+1) = minR;                % Uppdate radius
        C.calculateGradient(rayInd);         % Uppdate shell gradient
        C.calculateGradient(rayInd+1);       % Uppdate next shell gradient
        end

        % Add analytics
        [C.analyticals{rayInd,:}] = deal(Ans, AnH, AnL, Atp, Abp, Aip, Aop, Aac, iter-1);
        
        % Uppdate source
        source.P(rayInd,:) = p1;
        source.V(rayInd,:) = v1;
        source.path{rayInd} = rayPath;
        source.phase(rayInd) = totalPhase;
        source.totalPath(rayInd) = totalPath;
    end % End of main trace function and main loop -----------------------%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    %                         Support functions                           %
    %                                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %---------------------------------------------------------------------%
                      % Integration step function %
                      
    function [p1,v1] = intStep(p0,v0,shellInd)
        r0 = norm(p0);                                % Distance from center
        gV = -p0./r0;                                 % Refractive index gradient-vector (Center direction-vector)
        gV = gV./norm(gV);                            % Normalization of center vector
        theta = acos(dot(v0,gV)./norm(v0));           % Angle between ray-vector and gradient-vector
        nr = C.getN(shellInd,p0);                     % Get refractive index in ray-point based on the current n-gradient
        r = nr/(sin(theta) * -C.gradient(shellInd));  % Calculate local radius of curvature
        if abs(dot(gV,v0)) == 1 || isnan(theta)       % Check if vectors are paralell 
            v1 = v0;
            p1 = p0 + ds.*v1;
        else
            a = sign(det([v0' gV']));       % Find if vectors are reght or left oriented
            psi = ds * a / r;               % Rotation angle
            R = @(psi) [cos(psi) -sin(psi); % Rotation matrix
                sin(psi) cos(psi)];
            v1 = R(psi)*v0';                % Rotate ray-vector
            v1 = v1';   
            or = p0' + r.*R(a*pi/2)*v0';    % Center of rotation
            p1 = R(psi)*(p0'-or) + or;      % Rotate ray-point around center of rotation
            p1 = p1';
        end
    end

    %---------------------------------------------------------------------%
                      % Single step phase shift %    
    function dP = phaseShift(p0,shellInd)       
        nr = C.getN(shellInd,(p0)); % Get refractive index in ray-point based on the current n-gradient
        dP = ds*(nr-C.n(1)); % Phase shift relative to n0 = C.n(1)
    end

    %---------------------------------------------------------------------%
                      % Back-project ray to center line %
    function rIp = projectRay(p,v)
        t = -p(2)/v(2);
        rIp = p(1) + t*v(1);
    end
end