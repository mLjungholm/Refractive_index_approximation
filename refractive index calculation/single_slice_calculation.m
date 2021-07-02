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
s = source2d_variable_spacing([sStart,ones(length(sStart),1).*r.*1.2],ones(length(sStart),2).*[0,-1],n0);

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
ds = 10^3;
[n,~,shellR] = findRefractiveIndex_rungeKutta(s,peakPos,dPhase(:,1),n0,ds);
% shellR = shellR(n>0);
% n = n(n>0);

% Compare with known trace
% fprintf(fileID,'%6s %12s\r\n','x','exp(x)');
formatSpec = 's1 total path: %.4f, stest total path: %.4f, Diff is :%4f %%\n';
fprintf(formatSpec,s.totalPath(1),testS.totalPath(1),(testS.totalPath(1)-s.totalPath(1))/s.totalPath(1)*100)
formatSpec = 's1 phase: %.4f, stest phase: %.4f, Diff is :%4f %%\n';
fprintf(formatSpec,s.phase(1),testS.phase(1),(testS.phase(1)-s.phase(1))/s.phase(1)*100)


% figure(4)
% hold on; title('Refractive Index'); xlabel('radius [m]'); ylabel('n')
% plot(shellR,n,'bx')
% plot(peakPos,nFunc(peakPos./r),'--k')
% plotLine([r, max(n)],[r,min(n)],'k',4)
% plot(peakPos,nf,'gx')










