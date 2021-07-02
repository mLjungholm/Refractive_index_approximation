% calculate the phase shift for any given ray tracing
% input: source, stop-line, width of source
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tracing the rays
% Source parameters
v = [0;-1];
p = [0;1.2];
nRays = 1000;
width = 2.2;
n0 = 1.3;
n1 = 1.45;
% r = 1;

% Initiate source
s = source2d(p', v', nRays, width, n0, 'half');

% Initiate ray trace
ds = 10^3;
tic
runge_kutta_trace_known_gradient(s,ds,'parabolic',n0,n1,1)
toc

% Plot result
% figure(1)
% hold on; axis equal; grid on
% title('Runge-Kutta ray-trace')
% plotCircle(1,2*pi,1)
% plotLine([-1.2 0],[1.2 0],'k',1)
% plotLine([0 -1.2],[0 1.2],'k',1)
% s.plotTrace(1,'r')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Calculating phaseShift
%
% Simulated parameters of object
lambda = 550*12^-9;
r = lambda*6*2;
xRange = [0,1];

% Calculating phaseshift values at center line
s.getBacktrace(0);
phaseShift = s.getPhaseShiftValues(lambda,r,xRange);
phaseShift = [phaseShift; 0 r];
[peakVal,peakPos,peaksN] = findPeaks(phaseShift(:,2),phaseShift(:,1),0.2);

% Visualization of phase shift
figure(2)
hold on; title('Phase shift')
xlabel('radius [m]')
plot(phaseShift(:,2),phaseShift(:,1))
scatter(peakPos,peakVal)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Testing phase shift profile gained from peaks
%               against the full trace

phaseDiff = s.getPhaseDiff(lambda,r);
peaksDiff = getPhaseShift(peakVal,lambda);
dPhase = [peaksDiff, peakPos];

figure(3)
hold on; grid on
plot(phaseDiff(:,1),phaseDiff(:,2),'b')
plot(peakPos,peaksDiff,'*')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a new source with the peak positions and trace. This is to be used
% as a comarison to the aproximation trace.

testS = source2d_variable_spacing([peakPos(2:end),ones(peaksN-1,1)*1.2*r],[0,-1],n0);

steps = 10^3;
runge_kutta_trace_known_gradient(testS,steps,'parabolic',n0,n1,r);
testS.getBacktrace(0);
% Plot result
figure(4)
hold on; axis equal; grid on
title('Runge-Kutta ray-trace')
plotCircle(r,2*pi,4)
plotLine([-1.2*r 0],[1.2*r 0],'k',4)
plotLine([0 -1.2*r],[0 1.2*r],'k',4)
testS.plotTrace(4,'r')
testS.plotBacktrack(4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function used to get the phase shift from the peaks
function peaksDiff = getPhaseShift(peakVal,lambda)
halfsteps = (0:1:(length(peakVal)-2))./2.*lambda;
if peakVal(end-1) > 0.9
    p = halfsteps(end) + (1-peakVal(end))*lambda/2;
elseif peakVal(end-1) < 0.1
    p = halfsteps(end) + peakVal(end)*lambda/2;
else
    % If the value is here there is some sort of error
    disp('Error in getting the phase shift profile')
end
peaksDiff = [halfsteps';p];
end


