function n = homogenious_approximation(source, m_phase_shift, m_phase_pos,n0)
% I could define the source in here and give the source as an output
% object.
nShells = source.nRays;
n = zeros(nShells+1,1);
n(1) = n0;
shellR = m_phase_pos;
% function nr = firstApproximation(rayInd)

[ip0,ip1,~] = circleIntersect(shellR(1),source.P(1,:),source.V(1,:));
dv = (ip1-ip0);
d = sqrt(dv(1)^2+dv(2)^2); % distance traveled in shell
n(2) = m_phase_shift(2)/d + n0;

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
    n(rayInd+1) = dPhase(end)/dShells(end) + n0;
end
end