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

if ~isnan(app.lens.mLines{lineInd}.leftPhaseMax)
    
app.MaxLayerEditField.Value = app.lens.mLines{lineInd}.leftPhaseMax;
app.MinLayerEditField.Value = app.lens.mLines{lineInd}.leftPhaseMin;
end
if ~isnan(app.lens.mLines{lineInd}.centerZone)
app.CenterExclusionEditField.Value = app.lens.mLines{lineInd}.centerZone;
end
if ~isnan(app.lens.mLines{lineInd}.edgeZone)
app.EdgeExclusionEditField.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.edgeZone;
end
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