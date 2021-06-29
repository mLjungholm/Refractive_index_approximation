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

% lambda = 550*10^-9;
% r = 3*10^-5; % [10um]
% xRange = [0,1];
% phaseShift = s.getPhaseShiftValues(lambda,r,xRange);
% [peakVal,peakPos,peaksN] = findPeaks(phaseShift(:,2),phaseShift(:,1),0.1);
% 
% 
% figure(1)
% hold on
% plot(phaseShift(:,2),phaseShift(:,1))
% scatter(peakPos,peakVal)
% 
% 
dPhase = ones(length(peakVal)-2,1).*((1:1:length(peakVal)-2).*0.5)';
dPhase = [0; dPhase; abs(peakVal(end-1)-peakVal(end)) + dPhase(end)];
dPhase = dPhase.*lambda;
% dPhase = [0; dPhase; (dPhase(end) + phaseShift(end,1))];
% figure(2)
% title('Measured Phase Shift')
% xlabel('radius [m]')
% ylabel('Wavelength half shift [a.u]')
% hold on
% plot(realpos(end-20:end),realphase(end-20:end),'k')
% plot(peakPos(1:2),dPhase(1:2))

sStart = peakPos;
sStart = sStart(2:end);
s2 = source2d_variable_spacing([sStart,ones(length(sStart),1).*r.*1.2],ones(length(sStart),2).*[0,-1],1);

% Plot ray starting position and shells
% figure(3); title('ray start'); xlabel('radius [m]');
% grid on; hold on; axis equal
% s2.plotP(3)
% for rInd = 1:length(sStart)
%     plotCircle(peakPos(rInd),3)
%     plotLine([-r 0],[r 0],'k',3)
%     plotLine([0 r],[0 -r],'k',3)
% end


% Create refractive index profile
% n0 = 1;
[n,~,shellR] = findRefractiveIndex(s2,peakPos,dPhase,nGradient(1),1000);
shellR = shellR(n>0);
n = n(n>0);

figure(4)
hold on; title('Refractive Index'); xlabel('radius [m]'); ylabel('n')
plot(shellR,n,'bx')
plot(peakPos,nFunc(peakPos./r),'--k')
plotLine([r, max(n)],[r,min(n)],'k',4)
% plot(peakPos,nf,'gx')










