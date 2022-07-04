function uppdateRefractiveIndexPlot(app)
lineInd = app.SamplingLineSpinner.Value;
if ~isnan(app.lens.mLines{lineInd}.lambda)
    app.WavelengthnmEditField_2.Value = app.lens.mLines{lineInd}.lambda;
end
if ~isnan(app.lens.mLines{lineInd}.n0)
    app.OutsiderefractiveindexEditField_2.Value = app.lens.mLines{lineInd}.n0;
end
app.UIAxes.NextPlot = 'add';
if ~isempty(app.lens.mLines{lineInd}.refractiveGradient)
    plot(app.UIAxes, app.lens.mLines{lineInd}.gradientD,...
        app.lens.mLines{lineInd}.refractiveGradient,'r')
end
app.UIAxes.XMinorGrid = 'on';
app.UIAxes.YMinorGrid = 'on';
end