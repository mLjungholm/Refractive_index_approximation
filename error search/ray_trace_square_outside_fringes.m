close all

plotRays = false;
lambda = 540*10^-9;
r = 20*10^-5;
v = [0;-1];
p = [0;1.2].*r;
nRays = 1*10^4;
% nRays = 10;
SourceHalfWidth = r*1.5;
n0 = 1.4;
n1 = 1.45;

% Initiate source----------------------------------------------------------
s = Source_2d(p', v', nRays, SourceHalfWidth);
sVis = Source_2d(p', v', 16, SourceHalfWidth*1.11);
P0 = s.p0;
V0 = s.v0;

% Initiate ray trace-------------------------------------------------------
tic
sTraced = traceSquare(P0,V0,n0,n1,r,lambda);
sTracedVis = traceSquare(sVis.p0,sVis.v0,n0,n1,r,lambda);
toc
s.path = sTraced.path;
s.P = sTraced.P;
s.projection = project2y0(sTraced.V,sTraced.P);

sVis.path = sTracedVis.path;
sVis.P = sTracedVis.P;
sVis.projection = project2y0(sTracedVis.V,sTracedVis.P);

amp = [(sin(sTraced.Phase - pi/2)+1)./2,s.projection(:,1)];
amp = sortrows(amp,2);
amp = smoothdata(amp,'movmedian',5);


%%

map = sTraced.centerDist ~= 0 & s.projection(:,1) > 0 & s.projection(:,1) < r;
relativePhaseShift = sTraced.centerDist(map)./lambda;
phasePos = amp(map,2);

map = sTraced.centerDist ~= 0 & s.projection(:,1) > 0 & s.projection(:,1) <r;
[peakVal,mPeakPos,nPeaks] = findPeaks(amp(map,2),amp(map,1),0.1);

mPhaseShift = getPhaseShiftfromPeaks(peakVal,mPeakPos,0);
 

%%
close all

figure(1)
set(gcf,'position',[400,50,1300,1000])
subplot(2,2,[2,4])
hold on; axis equal; grid minor
sVis.plotRays('r')
sVis.plotProjection('--b')
plotCircle(r,2*pi)
plotLine([-1.2*r 0],[1.2*r 0],'k')
plotLine([0 -1.2*r],[0 1.2*r],'k')
title('Ray-Trace (homogenious refracitve index)')
legend('Ray Path','Ray Back Projection')

subplot(2,2,1)
hold on; grid minor
map = amp(:,2)>0 & amp(:,2)<r; 
plot(amp(map,2),amp(map,1),'r')
plot([0;0],[1;0],'k')
plot([r;r],[1;0],'k')
scatter(mPeakPos,peakVal,'bo')
title('Generated Interference Pattern')
ylabel('Interference amplitude [A.U]')
xlabel('radial distance [m]')

subplot(2,2,3)
hold on; grid minor
plot(phasePos,relativePhaseShift,'r')
plot(mPeakPos,mPhaseShift,'bo')
plot(mPeakPos,mPhaseShift,'b')
plot([r;r],[max(relativePhaseShift);0],'k')
legend('True phase Shift','Measured phase shift')
ylabel('$\Delta\varphi / \lambda$','Interpreter','latex')
xlabel('radial distance [m]')
title('True and Measured PhaseShift')


%%
% Ray trace function-------------------------------------------------------
function sTraced = traceSquare(P0,V0,n0,n1,r,lambda)
rayNums = size(P0,1);
sourcePath = cell(rayNums,1);
totalPath = zeros(rayNums,1);
phase = zeros(rayNums,1);
V = zeros(rayNums,2);
P = zeros(rayNums,2);
centerPoint = zeros(rayNums,2);
centerDist = zeros(rayNums,1);
for i = 1:rayNums
    rayPath = zeros(4,2).*nan;
    rayPath(1,:) = P0(i,:);
    [p0,p02,intersect] = circleIntersect(r,P0(i,:),V0(i,:));
    if ~intersect || isequal(p0, p02)
        sourcePath{i} = [P0(i,:); P0(i,1) -r];
        totalPath(i) = 0;
        phase(i) = 0;
        V(i,:) = V0(i,:);
        P(i,:) = [P0(i,1) -r];
        continue
    end
    rayPath(2,:) = p0;
    v0 = V0(i,:);
    N = p0./norm(p0);
    [v1, ~] = Snell(v0, N, n0, n1);
    cd = -p0(2)/v1(2);
    centerPoint(i,:) = [p0(1) + cd*v1(1), 0];
    centerDist(i) =  norm(rayPath(2,:)-rayPath(1,:))*n0...
        + norm(rayPath(2,:)-centerPoint(i,:))*n1 -norm(P0(1,2)).*n0;
    [~,p1,~] = circleIntersect(r,p0,v1);
    rayPath(3,:) = p1;
    N = p1./norm(p1);
    [v1, ~] = Snell(v1, N, n1, n0);
    d = -(r+p1(2))/v1(2);
    rayPath(4,:) = [p1(1) + d*v1(1), -r];
    totalPath(i) = norm(rayPath(2,:)-rayPath(1,:))*n0...
        + norm(rayPath(3,:)-rayPath(2,:))*n1...
        + norm(rayPath(4,:)-rayPath(3,:))*n1;
    phase(i) = centerDist(i)*2*pi/lambda;
    sourcePath{i} = rayPath;
    V(i,:) = v1;
    P(i,:) = rayPath(4,:);
end
sTraced.rayNums = rayNums;
sTraced.path = sourcePath;
sTraced.totalPath = totalPath;
sTraced.Phase = phase;
sTraced.V = V;
sTraced.P = P;
sTraced.P0 = P0;
sTraced.centerDist = centerDist;
sTraced.centerPoint = centerPoint;
end

function p1 = project2y0(v,p0)
rayNums = size(v,1);
p1 = zeros(rayNums,2);
for rayInd = 1:rayNums
    d = -p0(rayInd,2)/v(rayInd,2);
    p1(rayInd,:) = [p0(rayInd,1) + d*v(rayInd,1), 0];
end
end

function phaseShift = getPhaseShiftfromPeaks(peakVal,peakPos,outsideVal)
flipped = false;
if peakPos(1) < peakPos(end)
    peakPos = flipud(peakPos);
    peakVal = flipud(peakVal);
    flipped = true;
end
peakNum = length(peakPos);
phaseShift = zeros(peakNum,1);
phaseShift(1) = 0;
if peakVal(2) < 0.9 && peakVal(2) > 0.1
    phaseShift = 'error';
    return
end
if outsideVal == 0 && peakVal(2) > 0.9
    phaseShift(2) = 1/2;
elseif outsideVal == 1 && peakVal(2) < 0.1
    phaseShift(2) = 1/2;
else
    phaseShift(2) = 1;
end
for pInd = 3:peakNum-1
    phaseShift(pInd) = 1/2*(pInd - 2) + phaseShift(2);
end
phaseShift(end) = abs(peakVal(end-1)-peakVal(end))*1/2 + phaseShift(end-1);
if flipped
    phaseShift = flipud(phaseShift);
end
end





