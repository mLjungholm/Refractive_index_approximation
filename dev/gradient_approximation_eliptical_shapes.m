clc
close all

%-------------------------------------------------------------------------%
                            % True gradient %
lambda = 550*10^-9;
r = lambda*55*2;
n0 = 1.3; 
n1 = 1.45;
k = (n0-n1)/r^2;
nFunc = @(x) k.*x.^2 + n1; %<- True n(r) (parabolic)
ra = r;
rb = r/4;

%-------------------------------------------------------------------------%

% Test source
sTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% Control source
% sLinear = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% Slice object
C = IndexSliceEllipse(mPhasePos,mPhaseShift,n0,rb/ra);

% Homogenious approximation
homogenious_approximation_elliptical(sTest,C); 

rayInd = 'all';
steps = 10^3;
print_analytical = false;
gradient_approximation_elliptical(sTest, C, steps, rayInd);

figure(1)
hold on; grid on
C.plotNCurve
C.plotN()
plot(linspace(r,0,100),nFunc(linspace(r,0,100)),'b')

figure(2)
hold on; grid on
C.plotPhase
plot(truePhasePos,truePhaseShift,'r')
% 
% analytical = C.analyticals;
% figure(3)
% hold on; grid on
% x = 1:length(analytical{1});
% nt = nFunc(ones(length(x),1).*C.r(2));
% plot(x,analytical{1,1},x,analytical{1,2},x,analytical{1,3},x,nt)
% legend("ns","nHigh","nLow","nTrue")
% 
% 
% figure(4)
% hold on; grid on
% x = 1:length(analytical{1});
% plot(x,analytical{rayInd,4},x,analytical{rayInd,5},x,analytical{rayInd,6},x,analytical{rayInd,7})
% legend("Total Phase","Back Phase","Inner Phase","Outer Phase")


% iterText = string(1:1:itterations(2));
% legendText = ["Parabolic","Linear",iterText];
% legend(legendText);
