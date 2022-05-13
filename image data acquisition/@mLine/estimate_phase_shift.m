% fix the function so that it creates half phase values for the starting
% layer.

function estimate_phase_shift(this)
this.PP(:) = 0;
% Order the start and end layer
startLayer = min([this.leftPhaseMin, this.leftPhaseMax]);
endLayer = max([this.leftPhaseMin, this.leftPhaseMax]);


span = startLayer:1:endLayer; % All included layers in the span
spanNums = length(span);    % Number of layers in the span
startPhase = linspace(0,1,spanNums);  % Starting phase for each layer (in fraction of wavelength)

if this.leftPhaseMax < this.leftPhaseMin
    startPhase = flip(startPhase);
end

figure(1)
hold on
grid minor
plot([0;0],[0;0.5],'k')
plot([0;0],[0;0],'k.')
% Loop through all layers in the span
for spanInd = 1:spanNums
    points_in_layer = this.L{span(spanInd)};
    pNums = length(points_in_layer);
    
    % if start = high & peak = high -> shift = 1-start
    % if start = high & peak = low -> shift = start
    % if start = low & peak = high -> shift = 1-start
    % if start = low & peak = low -> shift = start
    p1Val = this.PV(points_in_layer(1)); % Amplitude of the first point
    if p1Val > 0.5
        startShift = (1-startPhase(spanInd))/2;
    else
        startShift = startPhase(spanInd)/2;
    end

    this.PP(points_in_layer(1)) = startShift;
    for pInd = 2:pNums
        cdelta = abs(this.PD(points_in_layer(pInd))-this.centerLine + this.leftEdge);
        if cdelta < this.centerZone
%             endVal = this.centerVal(spanInd);
%             peVal = this.PV(points_in_layer(end));
%             if peVal > 0.5
%                 phaseTemp = (pInd-2)/2 + startShift + (1-endVal)/2;
%             else
%                 phaseTemp = (pInd-2)/2 + startShift + endVal/2;
%             end
%             this.PP(points_in_layer(pInd)) = phaseTemp;
%             txt = num2str(span(spanInd));
%             d = this.PD(points_in_layer(end));
%             text(d,phaseTemp,txt)
%             plot([0;0],[0;0],'k.')
%             plot(d,phaseTemp,'b.','markerSize',0.1)
            continue
        else
            phaseTemp = (pInd-1)/2 + startShift;
            this.PP(points_in_layer(pInd)) = phaseTemp;
        end
    end
end
end
