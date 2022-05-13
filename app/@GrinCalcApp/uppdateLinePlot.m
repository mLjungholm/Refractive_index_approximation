function uppdateLinePlot(app)
app.UIAxesLine.NextPlot = 'replaceChildren';
plot(app.UIAxesLine,0,0)

if isempty(app.lens.mLines)
    return
end

if app.RawButton.Value
    interpolation = 'raw';
elseif app.GaussianButton.Value
    interpolation = 'gauss';
else
    interpolation = 'sgolay';
end

if app.ShowalllayersCheckBox.Value
    layerIndex = 'all';
    app.UIAxesLine.NextPlot = 'add';
else
    layerIndex = app.LayerSpinnerLine.Value;
    app.UIAxesLine.NextPlot = 'replaceChildren';
end

currentLineId = app.SamplingLineSpinner.Value;
app.lens.mLines{currentLineId}.AppPlotData(app.UIAxesLine,interpolation,layerIndex)
app.UIAxesLine.XMinorGrid = 'on';
app.UIAxesLine.YMinorGrid = 'on';



end