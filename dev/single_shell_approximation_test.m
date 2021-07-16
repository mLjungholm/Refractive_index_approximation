% Script for testing the refractive index aproximation of a single shell.
% The approximation should end upp with a linear gradient simmilar to the
% liniarization of the true parabolic function given by "nFunc"


% Testcase:
% n0 = 1.3
% n1 = 1.45
% nFunc = @(x) k.*x.^2 + n1;
% r = 12*lambda, lambda = 550*10^-9 


k = (n0-n1)/r^2;
shellInd = 1;
nApprox = @(x) shellGradient(shellInd,1).*(x.*r-shellR(shellInd)) + shellGradient(shellInd,2);
nFunc = @(x) k.*(x.*r).^2 + n1;
rx = linspace(shellR(shellInd+1),shellR(shellInd),100)./r;
% rx = linspace(0,shellR(shellInd),100)./r;
na = nFunc(shellR(shellInd)/r);
nb = nFunc(shellR(shellInd+1)/r);
ra = shellR(shellInd);
rb = shellR(shellInd+1);
kL = (na-nb)/(ra-rb);
nFuncL = @(x) kL.*(x.*r-ra) + na;

close all

figure(1)
titleText = strcat('Gradients for shell :',num2str(shellInd));
title(titleText)
hold on; grid on
plot(rx,nFunc(rx),'b')
plot(rx,nApprox(rx),'r')
plot(rx,nFuncL(rx),'g')
legend('True','Approximation','True-Linear')

%%
% now testing that a linear aproximation of the true parabolic refractive
% index profile results in a ray path and phase similar to the true
% example.

sTrue = Source_2d([peakPos(2),ones(1,1)*1.2*r],[0,-1]);
sLinear = Source_2d([peakPos(2),ones(1,1)*1.2*r],[0,-1]);
ds = 10^3;
ray_trace(sTrue,ds,n0,n1,r,'parabolic','meggit')
sTrue.projectRays([0,0],[1,0],'back');
ray_trace(sLinear,ds,n0,nFuncL(0),r,'linear','meggit')
sLinear.projectRays([0,0],[1,0],'back');

figure(2)
hold on; axis equal; grid on
title('Traced rays for the true,linear and aproximated gradient')
plotCircle(r,2*pi,2)
plotLine([-1.2*r 0],[1.2*r 0],'k')
plotLine([0 -1.2*r],[0 1.2*r],'k')
sTrue.plotRays('r')
sTrue.plotProjection('--r');
sLinear.plotRays('g')
sLinear.plotProjection('--g');

%%
% Now tracing the testsource using the alogrithm in the aproximation script

sSliceTest = Source_2d([peakPos(2),ones(1,1)*1.2*r],[0,-1]);
slice_trace_test(sSliceTest,ds,r,n0,kL,1)
sSliceTest.projectRays([0,0],[1,0],'back')

sSliceTest.plotRays('b')
sSliceTest.plotProjection('--b');
% legend('True','Approximation','True-Linear')

fprintf('\n\n\nEnd values for the tested rays \n')
scaleVal = 10^-7;
sTrue.getEndVals_single(1,'print',scaleVal);
sLinear.getEndVals_single(1,'print',scaleVal);
[phaseT,totalpathT,vT,pT] = sSliceTest.getEndVals_single(1,'print',scaleVal);

fprintf('\nsTrue n avrg  ->     n = %.4f \n',sTrue.phase(1)/sTrue.totalPath(1) + n0);
fprintf('sLinear n avrg ->    n = %.4f \n',sLinear.phase(1)/sLinear.totalPath(1) + n0);
fprintf('sSliceTest n avrg -> n = %.4f \n',sSliceTest.phase(1)/sSliceTest.totalPath(1) + n0);
fprintf('Max n-index for shell 1 is nMax = %.4f \n',nb);
fprintf('Avr n-index for shell 1 is nAvr = %.4f \n',(na+nb)/2);

% Using the getLinearGradient function we would get a max n of
ns = sSliceTest.phase(1)/sSliceTest.totalPath(1) + n0;
nmax = getLinearGradient_CurvedTrace(ra,rb,n0,ns,sSliceTest.phase(1),sSliceTest.phase(1));
fprintf('\nMax n-index using getLinearGradient is nMax = %.4f \n\n',nmax);
% this turn sout to be useless since it returns the same value as the
% average using the total phase. But setting the average value as the max n
% would result in a to low gradient.

% Using the  find total Phase function. I.e the back projection to the
% center of the lens and taking the measured phase value the test ray would
% get a phase of.

[phaseTot,~] = findTotPhase(pT,vT,peakPos,mPhaseShift);
fprintf('The phase from the back projection would be: phase = %.4f \n',phaseTot/scaleVal)
fprintf('Which would result in a n-Average of n = %.4f \n\n',phaseTot/totalpathT + n0);

figure(3)
grid on; hold on
title('Phase shift profiles')
plot(peakPos,mPhaseShift,'b')
plot(truePhasePos,truePhaseShift,'r'),
legend('True','Measured')

% Now try doing an itteration and check the values of both the traced phase
% and the measured projected phase.

function [phaseTot,xip] = findTotPhase(p0,v0,peakPos,peaksDiff)
t = -p0(2)/v0(2);
xip = p0(1) + t*v0(1); % x-axis intersection distance
% find colsest phaseShift values and interpolate
pind = find((peakPos-xip)>0,1,'last');
if isempty(pind)
    pind = 1;
elseif pind == length(peakPos)
    pind = length(peakPos) - 1;
end
p0 = peakPos(pind);
p1 = peakPos(pind+1);
shift0 = peaksDiff(pind);
shift1 = peaksDiff(pind+1);

phaseTot = (shift1-shift0)/(p1-p0)*(xip-p0) + shift0;
end
