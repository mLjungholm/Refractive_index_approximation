function uppdateLinePlot(app)

app.UIAxes.NextPlot = 'replaceChildren';
cla(app.UIAxes,'reset') 
% plot(app.UIAxes,0,0)
app.UIAxes.NextPlot = 'add';

if isempty(app.lens.mLines)
    return
end

if app.RawButton.Value
    interpolation = 'raw';
elseif app.GaussianButton.Value
    interpolation = 'gauss';
elseif app.SavGoButton.Value
    interpolation = 'sgolay';
else
    interpolation = 'all';
end

if app.ShowalllayersCheckBox.Value
    layerIndex = 'all';
    
else
    layerIndex = app.LayerSpinner.Value;
end

currentLineId = app.SamplingLineSpinner.Value;
app.lens.mLines{currentLineId}.AppPlotData(app.UIAxes,interpolation,layerIndex)
app.UIAxes.XMinorGrid = 'on';
app.UIAxes.YMinorGrid = 'on';

if app.ShowPeaks.Value
    app.lens.mLines{currentLineId}.AppPlotPeaks(app.UIAxes,interpolation,layerIndex)
end

end