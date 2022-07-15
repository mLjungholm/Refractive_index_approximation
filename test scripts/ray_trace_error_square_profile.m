% calculate the phase shift for any given ray tracing
% input: source, stop-line, width of source
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tracing the rays
% Source parameters
lambda = 545*10^-9;
rn = lambda*55*2;
v = [0;-1];
p = [0;1.2];
nRays = 1000;
width = rn;
n0 = 1.4072;
n1 = 1.473;
ra = rn;
rb = ra/4;
nProfile = 'parabolic';
% k = (n0-n1)/r^2;
% nFunc = @(x) k.*x.^2 + n1;

% k = (n0-n1)/rn^2;
nFunc = @(r,r0) (n0-n1)/rn^2.*r.^2 + n1;

% Initiate source
s = Source_2d(p', v', nRays, width);
sIm = Source_2d(p', v', 10, width);
% s2 =Source_2d(p', v', nRays, width);
% Initiate ray trace

ds = 10^3;
tic
ray_trace_elliptical(s,ds,n0,n1,ra,rb,nProfile,'meggit')
ray_trace_elliptical(sIm,ds,n0,n1,ra,rb,nProfile,'meggit')
s.projectRays([0,0],[1,0],'back');
toc


% Square trace-------------------------------------------------------------
% Initiate source
% s = Source_2d(p', v', nRays, width);
% sIm = Source_2d(p', v', 10, width);
% 
% % Initiate ray trace
% ds = 10^3;
% tic
% ray_trace(s,ds,n0,n1,r,nProfile,'snell')
% ray_trace(sIm,ds,n0,n1,r,nProfile,'snell')
% s.projectRays([0,0],[1,0],'back');
% toc
%--------------------------------------------------------------------------

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
n0 = 1.4072;

% Test source
sTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*rn],[0 -1]);

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
plot(linspace(rn,0,100),nFunc(linspace(rn,0,100),rn),'b')


%%

rn = C.r(1);
v = [0;-1];
p = [0;1.2];
nRays = 1000;
width = rn;
n0 = C.n(1);

% Initiate source
sVal = Source_2d(p', v', 10, width);
sVal2 = Source_2d(p', v', nRays, width);
gridNums = 10^3;
grin = GRIN2d_rotSym(C.n, C.r, gridNums);
% Initiate ray trace
stepSize = rn/10^3;
rayTrace2dGRIN_parallel(sVal,grin,stepSize,rn,n0)
rayTrace2dGRIN_parallel(sVal2,grin,stepSize,rn,n0)


%%
% Plot result
figure(1)
hold on; axis equal; grid minor
% % title('Runge-Kutta ray-trace')
plotCircle(rn,2*pi)
plotLine([-1.2*rn 0],[1.2*rn 0],'k')
plotLine([0 -1.2*rn],[0 1.2*rn],'k')
sVal.plotRays('r')
sIm.plotRays('b')

sVal2.projectRays([0,0],[1,0],'back');
[truePhaseShift, truePhasePos, relativePhaseShift] = create_interference_pattern(sVal2,lambda,nan,0);
figure(2)
% hold on; grid on
plot(truePhasePos,relativePhaseShift)
