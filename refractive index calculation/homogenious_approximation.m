% In use
function n = homogenious_approximation(source, C)
% I could define the source in here and give the source as an output
% object.
nShells = source.nRays;
n = zeros(nShells+1,1);
n(1) = C.n0;
m_phase_shift = C.mPhaseShift;
shellR = C.mPhaseRadius;
% function nr = firstApproximation(rayInd)

[ip0,ip1,~] = circleIntersect(shellR(1),source.P(1,:),source.V(1,:));
dv = (ip1-ip0);
d = sqrt(dv(1)^2+dv(2)^2); % distance traveled in shell
n(2) = m_phase_shift(2)/d + n(1);

for rayInd = 2:source.nRays
    dShells = zeros(rayInd,1); % Distance traveled in each shell
    for shellInd = 1:rayInd
        [ip0,ip1,~] = circleIntersect(shellR(shellInd),source.P(rayInd,:),source.V(rayInd,:));
        dv = (ip1-ip0);
        d = sqrt(dv(1)^2+dv(2)^2); % distance traveled in shell
        dShells(shellInd) = d;
    end
    % Substract the distance from inner shells.
    for shellInd = 1:rayInd-1
        dShells(shellInd) = dShells(shellInd) - dShells(shellInd+1);
    end
    
    dPhase = zeros(rayInd-1,1);
    
    for shellInd = 1:rayInd-1
        dPhase(shellInd) = dShells(shellInd)*(n(1+shellInd) - n(1));
    end
    
    
    % % dPhase = m_phase_shift(rayInd+1) - sum(dPhase);
    dPhase(rayInd) = m_phase_shift(rayInd+1) - sum(dPhase);
    n(rayInd+1) = dPhase(end)/dShells(end) + n(1);
end
C.n = n;
C.nh = n;
end