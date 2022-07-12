% calculate the phase shift for any given ray tracing
% input: source, stop-line, width of source
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tracing the rays
% Source parameters
lambda = 545*10^-9;
r = lambda*55*2;
v = [0;-1];
p = [0;1.2];
nRays = 1000;
width = r;
n0 = 1.4072;
n1 = 1.473;
ra = r;
rb = ra/4;
nProfile = 'square';

% Initiate source
s = Source_2d(p', v', nRays, width);

% Initiate ray trace
ds = 10^3;
tic
ray_trace(s,ds,n0,n1,r,nProfile,'snell')
s.projectRays([0,0],[1,0],'back');
toc

% Plot result
% figure(1)
% hold on; axis equal; grid minor
% % % title('Runge-Kutta ray-trace')
% plotCircle(r,2*pi)
% plotLine([-1.2*r 0],[1.2*r 0],'k')
% plotLine([0 -1.2*r],[0 1.2*r],'k')
% s.plotRays('r')
% s.plotEndRays(r/5)
% s.plotProjection('--b')


% s.getEndVals_single(9,'print');
% s2.getEndVals_single(9,'print');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Calculating phaseShift
%
% Simulated parameters of object

% r = 1;
% lambda = r/60;
% xRange = [0,1];
[truePhaseShift, truePhasePos, relativePhaseShift] = create_interference_pattern(s,lambda,nan,0);
% 
figure(2)
hold on; grid on
plot(truePhasePos,relativePhaseShift)
% 
[peakVal,mPhasePos,nPeaks] = findPeaks(truePhasePos,relativePhaseShift,0.2);
% peakPos(1) = r;
scatter(mPhasePos,peakVal,'bo')
% 
mPhaseShift = get_phase_shift_from_peaks(peakVal,lambda);
% 
figure(3)
hold on; grid on
plot(truePhasePos,truePhaseShift,'b')
plot(mPhasePos,mPhaseShift,'r')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lambda = 545*10^-9;
n0 = 1.37;

% Test source
sTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

C = IndexSliceEllipse(mPhasePos,mPhaseShift,n0,1);

% Homogenious approximation
homogenious_approximation_elliptical(sTest,C); 

rayInd = 'all';
steps = 10^3;
print_analytical = false;
gradient_approximation_elliptical(sTest, C, steps, rayInd);

figure(4)
hold on; grid on
C.plotNCurve
C.plotN()
% plot(linspace(r,0,100),nFunc(linspace(r,0,100)),'b')


%%

