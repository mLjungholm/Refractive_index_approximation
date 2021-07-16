% calculate the phase shift for any given ray tracing
% input: source, stop-line, width of source
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tracing the rays
% Source parameters
v = [0;-1];
p = [0;1.2];
nRays = 1000;
width = 1.5;
n0 = 1.3;
n1 = 1.45;
r = 1;
nProfile = 'parabolic';

% Initiate source
s = Source_2d(p', v', nRays, width);
% s2 =Source_2d(p', v', nRays, width);
% Initiate ray trace

ds = 10^3;
tic
ray_trace(s,ds,n0,n1,r,nProfile,'meggit')
s.projectRays([0,0],[1,0],'back');
toc
% tic
% ray_trace(s2,ds,n0,n1,r,nProfile,'rk')
% toc

% Plot result
% figure(1)
% hold on; axis equal; grid on
% title('Runge-Kutta ray-trace')
% plotCircle(1,2*pi,1)
% plotLine([-1.2 0],[1.2 0],'k')
% plotLine([0 -1.2],[0 1.2],'k')
% s.plotRays('r')
% s.plotProjection('--b')
% s2.plotRays('b')

% s.getEndVals_single(9,'print');
% s2.getEndVals_single(9,'print');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Calculating phaseShift
%
% Simulated parameters of object
lambda = 550*10^-9;
r = lambda*12*2;
xRange = [0,1];
[truePhaseShift, truePhasePos, relativePhasePos] = create_interference_pattern(s,lambda,r,1);

figure(2)
hold on; grid on
plot(truePhasePos,relativePhasePos)

[peakVal,peakPos,nPeaks] = findPeaks(truePhasePos,relativePhasePos,0.2);
scatter(peakPos,peakVal,'bo')

mPhaseShift = get_phase_shift_from_peaks(peakVal,lambda);

figure(3)
hold on; grid on
plot(truePhasePos,truePhaseShift,'b')
plot(peakPos,mPhaseShift,'r')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a new source with the peak positions and trace. This is to be used
% as a comarison to the aproximation trace.

sContr = Source_2d([peakPos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

ds = 10^3;
ray_trace(sContr,ds,n0,n1,r,nProfile,'meggit')
sContr.projectRays([0,0],[1,0],'back');

% Plot result
figure(4)
hold on; axis equal; grid on
plotCircle(r,2*pi,4)
plotLine([-1.2*r 0],[1.2*r 0],'k')
plotLine([0 -1.2*r],[0 1.2*r],'k')
sContr.plotRays('r')
sContr.plotProjection('--b')















