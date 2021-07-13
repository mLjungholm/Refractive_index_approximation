% Development script for the gradient aproximation algorithm

% Load data or generate new interference data with the 
% "createPhaseShift.m" script.

% create a source with variable spacing to match the positions of the
% measured phase values.

startP = [1/2 1.2];
startV = [0 -1];
nRays = 1;
width = 1;
testS = Source_2d(startP,startV,nRays,width);

lambda = 550*10^-9;
r = 12*lambda;

n = approximate_refracitve_index_first(source, phase_val, phase_pos,lambda,r);


function n = approximate_refracitve_index_first(source, phase_val, phase_pos)
% I could define the source in here and give the source as an output
% object.

nShells = source.nRays;
shellR = phase_pos;
n = zeros(nShells,1);

% function nr = firstApproximation(rayInd)
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
dPhase(rayInd) = ps(rayInd+1) - sum(dPhase);
nr = dPhase(end)/dShells(end) + n0;
end