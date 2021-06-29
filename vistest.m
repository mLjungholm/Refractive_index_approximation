realphase = s.phase;
realpos = s.backTrace(:,1);
realpos = realpos(realphase > 0);
realphase = realphase(realphase > 0);
realphase = realphase(realpos >= 0);
realpos = realpos(realpos >= 0);

realphase = realphase.*r;
realpos = realpos.*r;

dPhase = ones(length(peakVal)-2,1).*((1:1:length(peakVal)-2).*0.5)';
dPhase = [0; dPhase; abs(peakVal(end-1)-peakVal(end)) + dPhase(end)];
dPhase = dPhase.*lambda;

close all

figure(1)
hold on
plot(realpos,realphase,'k')
plot(peakPos,dPhase,'b')

