function uppdateLinePlot(app)

app.UIAxesLine.NextPlot = 'replaceChildren';
plot(app.UIAxesLine,0,0)
app.UIAxesLine.NextPlot = 'add';

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
    layerIndex = app.LayerSpinnerLine.Value;
end

currentLineId = app.SamplingLineSpinner.Value;
app.lens.mLines{currentLineId}.AppPlotData(app.UIAxesLine,interpolation,layerIndex)
app.UIAxesLine.XMinorGrid = 'on';
app.UIAxesLine.YMinorGrid = 'on';

if app.ShowPeaks.Value
    app.lens.mLines{currentLineId}.AppPlotPeaks(app.UIAxesLine,interpolation,layerIndex)
end

end