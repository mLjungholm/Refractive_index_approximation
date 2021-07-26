lc
close all

%-------------------------------------------------------------------------%
                            % True gradient %
lambda = 550*10^-9;
r = lambda*30*2;
n0 = 1.3; 
n1 = 1.45;
% k = (n0-n1)/r^2;
% nFunc = @(x) k.*x.^2 + n1; %<- True n(r) (parabolic)
% k = (n0-n1)/r;
% nFunc = @(x) k.*x + n1; %<- True n(r) (parabolic)
nFunc = @(x) n1*(1 - sign(x-1))/2;

%-------------------------------------------------------------------------%
                    % Refractive index approximations %
% Test source
sTest = Source_2d([mPhasePos(2:end),ones(nPeaks-1,1)*1.2*r],[0 -1]);

% Slice object
C = IndexSlice(mPhasePos,mPhaseShift,n0);

rayInd = 'all';
steps = 10^3;
print_analytical = false;

% Homogenious approximation
homogenious_approximation(sTest,C);

% Gradient approximation
gradient_approximation(sTest, C, steps, rayInd);

%-------------------------------------------------------------------------%
                         % Plotting values %
figure(1)
hold on; grid on
C.plotNCurve
C.plotN()
plot(linspace(r,0,100),nFunc(linspace(r,0,100)),'b')

% figure(3)
% hold on; grid on
% x = 1:length(analytical{1});
% plot(x,analytical{1},x,analytical{2},x,analytical{3},x,analytical{4},x,analytical{5})
% legend("ns","nHigh","nLow","nTrue","n0")
% 
% figure(4)
% hold on; grid on
% x = 1:length(analytical{1});
% plot(x,analytical{6},x,analytical{7},x,analytical{8},x,analytical{9})
% legend("Total Phase","Back Phase","Inner Phase","Outer Phase")