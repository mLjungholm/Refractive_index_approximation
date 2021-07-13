% Refractive index calculation for one slice
% In this version the slice needs to be circularly cymetric

% Input:    phase shift values
%           wavelength
%           slice dimensions
%           medium refractive index


% Import phase shift values
% Names:    phasePos - phase shift position coordinates
%           phaseShift - full integer phase shift values for positions

% Start with outer most ray -> check the orientation of the positions and
% phase shift values

% Define outer shell
close all

% figure(1)
% hold on
% plot(phaseShift(:,2),phaseShift(:,1))
% scatter(peakPos,peakVal)


sStart = peakPos;
sStart = sStart(2:end);
sTest = source2d_variable_spacing([sStart,ones(length(sStart),1).*r.*1.2],ones(length(sStart),2).*[0,-1],n0);

% Plot ray starting position and shells
% figure(3); title('ray start'); xlabel('radius [m]');
% grid on; hold on; axis equal
% s2.plotP(3)
% for rInd = 1:length(sStart)
%     plotCircle(peakPos(rInd),3)
%     plotLine([-r 0],[r 0],'k',3)
%     plotLine([0 r],[0 -r],'k',3)
% end
k = (n0-n1)/r^2;
nFunc = @(r) k.*r.^2 + n1;

% Create refractive index profile
n0 = 1.3;
ds = 10^3;
[n,shellGradient,shellR] = findRefractiveIndex(sTest,peakPos,dPhase(:,1),n0,ds);
% shellR = shellR(n>0);
% n = n(n>0);

% Compare with known trace
% fprintf(fileID,'%6s %12s\r\n','x','exp(x)');
% formatSpec = '\n \nsTest total path: %.4f \nsContr total path: %.4f Diff is :%4f %%\n';
% fprintf(formatSpec,sTest.totalPath(1)/r,sContr.totalPath(1)/r,(sContr.totalPath(1)-sTest.totalPath(1))/sTest.totalPath(1)*100)
% formatSpec = 'sTest phase: %.4f \nsContr phase: %.4f Diff is :%4f %%\n';
% fprintf(formatSpec,sTest.phase(1)*10^7,sContr.phase(1)*10^7,(sContr.phase(1)-sTest.phase(1))/sTest.phase(1)*100)
% formatSpec = 'sTest average n: %.4f \nsTest (phase/d + n0) = %.4f \nsContr average n: %.4f  \n';
% fprintf(formatSpec,n(2), sTest.phase(1)/sTest.totalPath(1) + n0,sContr.phase(1)/sContr.totalPath(1) + n0)


% figure(4)
% hold on; title('Refractive Index'); xlabel('radius [m]'); ylabel('n')
% plot(shellR,n,'bx')
% plot(peakPos,nFunc(peakPos),'--k')
% plotLine([r, max(n)],[r,min(n)],'k',4)
% plot(peakPos,n,'gx')

% figure(1)
% % yyaxis left
% sContr.plotRay(1,'-r',1)
% xi = sContr.plotRay_back(1,0,'--r',1);
% xii = sTest.plotRay_back(1,0,'--b',1);

% figure(2)
% hold on; grid on
% plot(phaseDiff(:,1),phaseDiff(:,2),'m')
% plot([xi;xi],[0;max(phaseDiff(:,2))],'-k')
% plot([xii;xii],[0;max(phaseDiff(:,2))],'-k')


