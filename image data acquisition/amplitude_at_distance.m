function amplitude_at_distance(d,testLine)
[~,closestInd] = min(abs(testLine.d - d));
amp = zeros(testLine.imNums,1);
for layerInd = 1:testLine.imNums
%     amp(layerInd) = sum(testLine.gaussPoints(span,layerInd))/3;
%     amp(layerInd) = testLine.gaussPoints(closestInd,layerInd)/max(testLine.gaussPoints(:,layerInd));
    amp(layerInd) = testLine.gaussPoints(closestInd,layerInd);
end

lmm = smoothdata(amp,'movmedian',4);
lg = smoothdata(amp,'gaussian',10);
ls = smooth(lmm,4,'sgolay');
x = 1:testLine.imNums;

% [psor,lsor] = findpeaks(lg,'MinPeakDistance',10,'SortStr','descend');

startL = testLine.gaussPoints(closestInd,1);
endL = testLine.gaussPoints(closestInd,testLine.imNums);
ampSpan = max(testLine.gaussPoints(closestInd,:))-min(testLine.gaussPoints(closestInd,:));
deltaVal = abs(round((startL-endL)/ampSpan,2));

figure(1)
subplot(2,1,1)
hold on
grid minor
plot(x,amp,'*','color',[0 0.4470 0.7410].*0.8)
plot(x,lg,'color','#FF0000')
plot(x,ls,'color','#EDB120')

plot([1;testLine.imNums],[1;1].*startL,'k');
plot([1;testLine.imNums],[1;1].*endL,'k');
plot([0.8;0.8].*testLine.imNums,[startL;endL],'r','linewidth',1.5)
labelText = strcat('\leftarrow \Delta = ',{' '},num2str(deltaVal));
labelPoint = startL-(startL-endL)/2;
text(0.8.*testLine.imNums,labelPoint,labelText)
text(testLine.imNums*1.01,startL*1.01,'start value');
text(testLine.imNums*1.01,endL*1.01,'end value');
title(strcat('Amplitude span of intensity values at distance:',{' '},num2str(d))) 
ylabel('Intensity [a.u]')
xlabel('Distance [pixels]')

subplot(2,1,2)
hold on
grid minor
i = 1;
plot(testLine.d,testLine.gaussPoints(:,i))
i = testLine.imNums;
plot(testLine.d,testLine.gaussPoints(:,i))
maxP = max(max(testLine.gaussPoints));
minP = min(min(testLine.gaussPoints));
plot(testLine.leftEdge.*[1;1],[maxP;minP],'k')
plot(testLine.rightEdge.*[1;1],[maxP;minP],'k')
plot(testLine.centerLine.*[1;1],[maxP;minP],'k')
plot(testLine.d(closestInd).*[1;1],[maxP;minP],'color','r')

end


