function clearOustidePeaks(app)
lineInd = app.SamplingLineSpinner.Value; 
layerMin = app.lens.mLines{lineInd}.leftPhaseMin;
layerMax = app.lens.mLines{lineInd}.leftPhaseMax;
startLayer = min([layerMin, layerMax]);
endLayer = max([layerMin, layerMax]);


% Remove peaks outside the span of the half shift layer
pinds = find(app.lens.mLines{lineInd}.Pinc);
for pi = 1:length(pinds)
    if app.lens.mLines{lineInd}.PL(pinds(pi)) < startLayer
        app.lens.mLines{lineInd}.Pinc(pinds(pi)) = false;
    elseif app.lens.mLines{lineInd}.PL(pinds(pi)) > endLayer
        app.lens.mLines{lineInd}.Pinc(pinds(pi)) = false;
    end
end

centerZone = app.lens.mLines{lineInd}.centerZone - app.lens.mLines{lineInd}.leftEdge;
edgeZone = app.lens.mLines{lineInd}.edgeZone - app.lens.mLines{lineInd}.leftEdge;
pinds = find(app.lens.mLines{lineInd}.Pinc);
if ~isnan(centerZone)
    for pi = 1:length(pinds)
        if app.lens.mLines{lineInd}.PD(pinds(pi)) > centerZone
            app.lens.mLines{lineInd}.Pinc(pinds(pi)) = false;
        end
    end
end
if ~isnan(edgeZone)
    for pi = 1:length(pinds)
        if app.lens.mLines{lineInd}.PD(pinds(pi)) < edgeZone
            app.lens.mLines{lineInd}.Pinc(pinds(pi)) = false;
        end
    end
end

end

% centerLine = nan;       % Coordinate for the center of the lens
%         centerLineIndex = nan;
%         leftEdge = nan;         % Edge of lens in relation to plot
%         leftEdgeIndex = nan;
%         rightEdge = nan;        % -||-
%         rightEdgeIndex = nan;
%         leftPhaseMax = nan;    % Index for the layer used as reference point for the phase shift
%         leftPhaseMin = nan;
%         
% %         sgolayZone = 0;
% %         centerVal = [];
%         
%         centerZone = nan;
% %         edgeZone = nan;
%         sgolayEdge = nan;
% %         sgolayCenter = nan;
%         
%         %         Variables from the trace peaks function
%         PD = [] % Point distance
%         PL = [] % Point Layer
%         PS = [] % Point series
%         PP = [] % Point Phase
%         PV = [] % peak value
%         L = {} % Each cell in L(layer) = [points];
%         S = {} % Series
%         Pinc = [] % Include point in calulations & plots
%         seriesNums = 0;