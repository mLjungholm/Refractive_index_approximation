% calculate the phase shift for any given ray tracing
% input: source, stop-line, width of source
close all
% 
% figure(1)
% scatter(s.backTrace(:,1),s.phase)

lambda = 550*10^-9;
% r = 1.65*10^-5; % [10um]
r = lambda*10*2;
xRange = [0,1];
% phaseShift = s.getPhaseShiftValues(lambda,r,xRange);

% phaseShift = s.phase.*scale;
% phasePos = s.backTrace(:,1);
% phaseShift = phaseShift(phasePos > -1);
% phasePos = phasePos(phasePos > -1);
% phaseShift = phaseShift(phasePos < 1);
% phasePos = phasePos(phasePos < 1);
% phasePos = phasePos.*scale;
% phaseFunc = @(x) (cos(x*2*pi/lambda + pi) + 1)/2;
% dphase = phaseFunc(phaseShift);
phaseShift = s.getPhaseShiftValues(lambda,r,xRange);
% f = fit(phaseShift(:,2),phaseShift(:,1),'smoothingspline','SmoothingParam',0.9);
[peakVal,peakPos,peaksN] = findPeaks(phaseShift(:,2),phaseShift(:,1),0.1);

figure(1)
hold on
plot(phaseShift(:,2),phaseShift(:,1))
scatter(peakPos,peakVal)



% 
% 
% function [peakVal2,peakPos2,nrPeaks] = findPeaks(coords,vals,threshold,minSpace)
% maxVal = max(vals);
% edgeP = [coords(1),coords(end)];
% edgeV = [vals(1), vals(end)];
% coords = coords(vals > maxVal*threshold);
% vals = vals(vals > maxVal*threshold);
% 
% peakVal = zeros(length(coords),1); peakVal(:) = nan;
% peakPos = zeros(length(coords),1); peakPos(:) = nan;
% y1 = vals(1);
% for ind = 2:length(coords)-1
%     y2 = vals(ind);
%     y3 = vals(ind + 1);
%     if y1 < y2 && y3 < y2
%         peakVal(ind) = y2;
%         peakPos(ind) = coords(ind);
%     elseif ind == length(coords)-1 && y3 > y2        
%         peakVal(ind+1) = y3;
%         peakPos(ind+1) = coords(ind+1);
%     end
%     y1 = y2;
% end
% 
% peakVal = peakVal(~isnan(peakVal));
% peakPos = peakPos(~isnan(peakPos));
% 
% peakVal2 = zeros(length(peakVal),1); peakVal2(:) = nan;
% peakPos2 = zeros(length(peakVal),1); peakPos2(:) = nan;
% endoflist = 0;
% startind = 1;
% peakind = 1;
% while ~endoflist
%     testList = zeros(length(peakVal),2).*nan;
%     x0 = peakPos(startind);
%     v0 = peakVal(startind);
%     testList(startind,:) = [x0,v0];
%     for ind = startind+1:length(peakVal)
%         x1 = peakPos(ind);
%         v = peakVal(ind);
%         d = abs(x1-x0);
%         if d < minSpace
%             testList(ind,:) = [x1,v];
%             x0 = x1;
%         else
%             startind = ind;
%             break
%         end
%     end
%     
%     [maxVal,indl] = max(testList(:,2));
%     peakVal2(peakind) = maxVal;
%     peakPos2(peakind) = testList(indl,1);
%     peakind = peakind + 1;
%     if ind == length(peakVal)
%         if startind == ind
%             x0 = peakPos(startind);
%             v0 = peakVal(startind);            
%             peakVal2(peakind) = v0;
%             peakPos2(peakind) = x0;
%         end
%         break
%     end
% end
% 
% 
% peakVal2 = peakVal2(~isnan(peakVal2));
% peakPos2 = peakPos2(~isnan(peakPos2));
% peakVal2 = [edgeV(1); peakVal2; edgeV(end)];
% peakPos2 = [edgeP(1); peakPos2; edgeP(end)];
% 
% nrPeaks = length(peakPos2);
% end

% %% This part is for visualization and creating the phase shift plots

% % Tried to make a function to find the peaks.
% % Did not work du to a bit crazy data. Have to do it manualy yyaaay....(sarcasm)
% % 
% function [peakVal,peakPos] = findPeaks(coords,vals)
% peakVal = zeros(length(coords),1); peakVal(:) = nan;
% peakPos = zeros(length(coords),1); peakPos(:) = nan;
% y1 = vals(1);
% for ind = 2:length(coords)-1
%     y2 = vals(ind);
%     y3 = vals(ind + 1);
%     if y1 < y2 && y3 < y2
%         peakVal(ind) = y2;
%         peakPos(ind) = coords(ind);
%     end
%     y1 = y2;
% end
% peakVal = peakVal(~isnan(peakVal));
% peakPos = peakPos(~isnan(peakPos));
% end



% xp = linspace(0,1,100);
% lambda = 5*10^-7;
% shift = xp.*lambda;
% phaseFunc = @(x) (cos(x*2*pi/lambda - pi) + 1)/2;
% dp = phaseFunc(shift);
% 
% figure(1)
% plot(xp,dp)