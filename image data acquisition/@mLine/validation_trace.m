function validation_trace(this)

% Set up variables
r0 = this.gradientD(1); % Maybe make a function that cheks in what order the values are and sorts them therafter.
v0 = [0 -1];
p0 = [0 1.2].*r0;
nRays = 1000;
raySteps = 1000;
gridNums = 1000;
width = r0;
n0 = this.n0;


% Create test Source
s = Source_2d(p0,v0, nRays, width);

% Create Grin system
grin = GRIN2d_rotSym(fliplr(this.refractiveGradient),fliplr(this.gradientD),gridNums);

% Trace rays
rayTrace2dGRIN_parallel(s,grin,r0/raySteps,r0,n0);
s.projectRays([0,0],[1,0],'back');

% Create interference pattern. 
[~, truePhasePos, relativePhaseShift] = create_interference_pattern(s,this.lambda,nan,0);

% [peakVal,mPhasePos,~] = findPeaks(truePhasePos,relativePhaseShift,0.2);

this.validS = s;
this.validP = truePhasePos;
this.validPS = relativePhaseShift;


% Run a trace with a different gradient


% Create a figure to compare phase shift for the estimated refractive index
% and the microscopy images.
figure(1)
hold on
layerIndex = this.leftPhaseMin;
grid minor
plot((max(this.validP)-this.validP)./this.pixelSize,this.validPS,'r')
plot(this.d,this.gaussPoints(:,layerIndex),'color',[0 0.4470 0.7410].*0.8)
plot(this.gaussPks{layerIndex}(:,2),this.gaussPks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.4)
this.plotSupportLines()


% mPhaseShift = get_phase_shift_from_peaks(peakVal,lambda);

end
