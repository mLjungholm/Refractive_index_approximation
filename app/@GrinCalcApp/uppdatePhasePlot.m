function uppdatePhasePlot(app)
app.UIAxes.NextPlot = 'replaceChildren';
cla(app.UIAxes,'reset')
app.UIAxes.NextPlot = 'add';
if app.PeakSelectButton.Value
    tools = 'peakselect';
elseif app.PeakTraceButton.Value    
    tools = 'peaktrace';
elseif app.PhaseGradientButton.Value
    tools = 'phasegradient';

end

layerIndex = app.LayerSpinner.Value;
lineInd = app.SamplingLineSpinner.Value;

if ~isnan(app.lens.mLines{lineInd}.layerPhaseStep)
    app.LayerPhaseShift.Value = app.lens.mLines{lineInd}.layerPhaseStep;
end
% Uppdate the peak table
if ~isempty(app.lens.mLines{lineInd}.S)
    pinds = find(app.lens.mLines{lineInd}.PS);
    group = app.lens.mLines{lineInd}.PS(pinds);
    layer = app.lens.mLines{lineInd}.PL(pinds);
    phase = app.lens.mLines{lineInd}.PP(pinds);
    pinc = app.lens.mLines{lineInd}.Pinc(pinds);
    tab = table(group,pinds,layer,phase,pinc);
    app.PhaseTable.Data = tab;
end

switch tools
    case 'peakselect'        
        app.lens.mLines{lineInd}.AppPlotData(app.UIAxes,'gauss',layerIndex)
        app.lens.mLines{lineInd}.AppPlotData(app.UIAxes,'sgolay',layerIndex)
%         app.lens.mLines{lineInd}.AppPlotPeaks(app.UIAxes,'gauss',layerIndex)
%         app.lens.mLines{lineInd}.AppPlotPeaks(app.UIAxes,'sgolay',layerIndex)
        app.lens.mLines{lineInd}.AppPlotPeaks(app.UIAxes,'pks',layerIndex)
        app.UIAxes.XMinorGrid = 'on';
        app.UIAxes.YMinorGrid = 'on';
        if ~isempty(app.lens.mLines{lineInd}.S)
%         if app.lens.mLines{lineInd}.seriesNums
%         if ~isempty(app.PhaseTable.Data)
            app.plotSelectedPeaks()
        end
    case 'peaktrace'
        if isempty(app.lens.mLines{lineInd}.S)
            app.addTwoConsole('Error: There are no selected peaks')
            app.PeakSelectButton.Value = true;
            return
        end
        app.plotSelectedPeaks()
        app.lens.mLines{lineInd}.AppPlotData(app.UIAxes,'gauss','all')
        app.lens.mLines{lineInd}.AppPlotPeaks(app.UIAxes,'pks','all')
        app.UIAxes.XMinorGrid = 'on';
        app.UIAxes.YMinorGrid = 'on';
    case 'phasegradient'
        if any(app.lens.mLines{lineInd}.PP)
            app.lens.mLines{lineInd}.AppPlotPhaseGradient(app.UIAxes)
            app.UIAxes.XMinorGrid = 'on';
            app.UIAxes.YMinorGrid = 'on';
        end
end

end
%
% if app.SamplingLineSpinnerPhase.Value == app.currentSamplingLine
%     return
% elseif app.SamplingLineSpinnerPhase.Value > length(app.lens.mLines)
%     app.SamplingLineSpinnerPhase.Value = app.currentSamplingLine;
%     return
% elseif isempty(app.lens.mLines{app.SamplingLineSpinnerPhase.Value})
%     app.SamplingLineSpinnerPhase.Value = app.currentSamplingLine;
%     return
% end
%
% app.currentSamplingLine = app.SamplingLineSpinnerPhase.Value;
% app.SamplingLineSpinner.Value = app.SamplingLineSpinnerPhase.Value;
%
% lineInd = app.SamplingLineSpinnerPhase.Value;
% if isempty(app.lens.mLines{lineInd}.PP)
%     return
% end
%
% for pInd = 1:app.lens.mLines{lineInd}.pksNums
%     if app.lens.mLines{lineInd}.PP(pInd) ~= 0
%         plot(app.UIAxesPhase,app.lens.mLines{lineInd}.PD(pInd),app.lens.mLines{lineInd}.PP(pInd),'b.','markerSize',15)
%     end
% end
% maxX = app.lens.mLines{lineInd}.centerLine;
% maxY = max(app.lens.mLines{lineInd}.PP);
% startLine = [0 0; 0 maxY];
% endLine = [maxX 0; maxX maxY];
% plot(app.UIAxesPhase, startLine(:,1),startLine(:,2),'k')
% plot(app.UIAxesPhase,endLine(:,1),endLine(:,2),'k')