function uppdateApp(app)
cla(app.UIAxes,'reset')
tabId = app.TabGroup.SelectedTab.Tag; 
switch tabId
    case 'stack'
        app.uppdateImstackImage()
    case 'line'
        app.uppdateLinePlot()
        app.uppdateLineData()
    case 'calibration'
        cla(app.UIAxes,'reset')
        if ~isempty(app.lens.calibIm)
            imshow(app.lens.calibIm,'Parent',app.UIAxes)
        end
    case 'phase'
        if isempty(app.LineDataTable.Data{3,2})...
                || isempty(app.LineDataTable.Data{1,2})...
                || isempty(app.LineDataTable.Data{2,2})...
                || isempty(app.lens.mLines{app.SamplingLineSpinner.Value}.gaussPks)
            app.addTwoConsole('The line data must be complete befor estimating the pahse')
            app.TabGroup.SelectedTab = app.TabGroup.Children(2);
            return
        end
        app.uppdatePhasePlot()
    case 'nindex'
        cla(app.UIAxes,'reset')
        if isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.pixelSize)
            app.addTwoConsole('The data must be calibrated before calculating the refractive index')
            app.TabGroup.SelectedTab = app.TabGroup.Children(5);
        elseif isempty(app.lens.mLines{app.SamplingLineSpinner.Value}.S)
            app.addTwoConsole('No peaks have been selected for estimation')
            app.TabGroup.SelectedTab = app.TabGroup.Children(2);
        elseif ~isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.refractiveGradient)
            app.lens.mLines{app.SamplingLineSpinner.Value}.AppPlotRefractiveGradient(app.UIAxes)
        else
            app.uppdateRefractiveIndexPlot()
        end
end
end


% >> h = app.TabGroup.Children(2);
% You can then directly assign it to SelectedTab property.
% >> app.TabGroup.SelectedTab = h;