function AppCalculateRefractiveIndex(this)
% Set up baseVals
r = this.edgeFit.*this.pixelSize;
% r = (this.centerLine-this.leftEdge).*this.pixelSize;

% Generate testPoints
% This should be equal to the number of measured phase points. This is not 
% needed since we are using the fitted profile.
pointNums = 20;
% dpoints = this.d - this.edgeFit;

nPeaks = pointNums;
mPhasePos = fliplr(linspace(0,this.edgeFit,pointNums))';

% for pind = 2:pointNums-1
%     closest = this.d-mPhasePos(pind)


% mPhasePos = fliplr(linspace(0,this.centerLine-this.leftEdge,pointNums))';
mPhaseShift = feval(this.phase_func,mPhasePos);
mPhasePos = mPhasePos.*this.pixelSize;
mPhaseShift(mPhaseShift < 0) = 0;
mPhaseShift = mPhaseShift.*this.lambda;
% mPhasePos = mPhasePos.*this.pixelSize;

% Create test source
sTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% Create slice object
C = IndexSliceEllipse(mPhasePos,mPhaseShift,this.n0,1);

% Homogenious approximation
homogenious_approximation_elliptical(sTest,C); 

rayInd = 'all';
steps = 10^3;
print_analytical = false;
gradient_approximation_elliptical(sTest, C, steps, rayInd);
this.testSlice = C;
this.refractiveGradient = C.n;
this.gradientD = C.r;

% figure(1)
% hold on; grid on
% C.plotNCurve
% C.plotN()
% xlabel('Position [m]')
% ylabel('Refractive index')
% plot(linspace(r,0,100),nFunc(linspace(r,0,100)),'b')

% figure(2)
% hold on; grid on
% C.plotPhase
end