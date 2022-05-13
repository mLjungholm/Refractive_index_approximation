% This function searches through the intensity values at the specified edge
% to find at what layer the edge is att a minimum or maximum. It returns
% the reference layers [min max], the total amplitude values of the curve at
% the reference edge and the difference between the start and end value of
% the span.

% Input:
% testLine - line to be tested
% edge - "left" or "right". What edge of the line to be tested
% plotPeaks - [0,1] Flag if to plot the results.

% Outputs:
% refLayer - [min max] index of the reference layers
% amplitude - [min max] min and max values of the span
% deltaVal - fractional difference between start and endpoint of the span

function [refLayer, amplitude, deltaVal] = find_reference_layer(testLine,edge,plotPeaks)

if isequal(edge,'left')
    [~,closestInd] = min(abs(testLine.d - testLine.leftEdge));
elseif isequal(edge,'right')
    [~,closestInd] = min(abs(testLine.d - testLine.rightEdge));
else
    fprintf('Error: Invalid reference edge \n')
    return
end

% span = ones(3,1).*closestInd + [-1;0;1];
edgeVal = zeros(testLine.imNums,1);

for layerInd = 1:testLine.imNums
    %     edgeVal(layerInd) = sum(testLine.gaussPoints(span,layerInd))/3;
    %     edgeVal(layerInd) = testLine.gaussPoints(closestInd,layerInd)/max(testLine.gaussPoints(:,layerInd));
    edgeVal(layerInd) = testLine.gaussPoints(closestInd,layerInd);
end

% lmm = smoothdata(edgeVal,'movmedian',4);
lg = smoothdata(edgeVal,'gaussian',10);

[~,peakPointMin] = findpeaks(-lg,'MinPeakDistance',10,'SortStr','descend');
[~,peakPointMax] = findpeaks(lg,'MinPeakDistance',10,'SortStr','descend');

refLayer = [peakPointMin(1) peakPointMax(1)];
amplitude = [min(edgeVal) max(edgeVal)];

startL = edgeVal(1);
endL = edgeVal(end);
deltaVal = abs(round((startL-endL)/(amplitude(1)-amplitude(2)),2));

if plotPeaks == 1
    x = 1:testLine.imNums;
    figure(1)
    hold on
    plot(x,edgeVal,'*','color',[0.9290 0.6940 0.1250].*0.8)
    plot(x,lg,'color',[0 0.4470 0.7410],'linewidth',1.5)
    plot(peakPointMin(1),lg(peakPointMin(1))-norm(amplitude)*0.02,'^','color', '#FF0000','MarkerFaceColor',[0.6350 0.0780 0.1840],'markersize',7)
    labelText = strcat('minimun at layer:',{' '},num2str(refLayer(1)));
    text(peakPointMin(1),lg(peakPointMin(1))-norm(amplitude)*0.055,labelText)
    plot(peakPointMax(1),lg(peakPointMax(1))+norm(amplitude)*0.02,'v','color', '#FF0000','MarkerFaceColor',[0.6350 0.0780 0.1840],'markersize',7)
    labelText = strcat('maximun at layer:',{' '},num2str(refLayer(2)));
    text(peakPointMax(1),lg(peakPointMax(1))+norm(amplitude)*0.055,labelText)
    plot([1;testLine.imNums],[1;1].*startL,'k');
    plot([1;testLine.imNums],[1;1].*endL,'k');
    plot([0.8;0.8].*testLine.imNums,[startL;endL],'color',[0.6350 0.0780 0.1840],'linewidth',1.2)
    labelText = strcat('\leftarrow \Delta = ',{' '},num2str(deltaVal));
    labelPoint = startL-(startL-endL)/2;
    text(0.8.*testLine.imNums,labelPoint,labelText)
    text(testLine.imNums*1.01,startL*1.01,'start value');
    text(testLine.imNums*1.01,endL*1.01,'end value');
    if isequal(edge,'left')
        titleString = 'Amplitude span of intensity values at the Left edge';
    elseif isequal(edge,'right')
        titleString = 'Amplitude span of intensity values at the Right edge';
    end
    title(titleString)
    ylabel('Intensity [a.u]')
    xlabel('Layer number')
    grid minor
end
end