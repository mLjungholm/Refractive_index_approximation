function checkImportedData(app)
if ~isequal(class(app.lens),'Imstack')
    app.addTwoConsole('Error: The imported data is not an GrinCalc project')
    app.lens = nan;
    return
end
if ~isempty(app.lens.imStack)
    app.LayerSpinner.Limits = [1,app.lens.imNums];
    app.LayerSpinnerLine.Limits = [1,app.lens.imNums];
    app.LayerSpinner.Enable = true;
    app.CropStackButton.Enable = true;
end
if ~isempty(app.lens.imStack)
    app.ShowCropedIm.Enable = true;
    app.showCropRec.Enable = true;
    app.showCropRec.Value = true;
    app.CreateSupportLineButton.Enable = true;
    app.CreateSamplingLineButton.Enable = true;
end
if ~isempty(app.lens.mLines)
    app.SamplingLinesCheckBox.Value = true;
    app.SamplingLinesCheckBox.Enable = true;
    app.DeleteLineButton.Enable = true;
    app.currentLineId = find(~cellfun(@isempty,app.lens.mLines),1);
    app.SamplingLineSpinner.Value = app.currentLineId;
    if length(app.lens.mLines) > 1
        app.SamplingLineSpinner.Limits = [1,length(app.lens.mLines)];
        app.SamplingLineSpinner.Enable = true;
    else
        app.SamplingLineSpinner.Value = app.currentLineId;
    end
    for ind = 1:length(app.lens.mLines)
        if ~isempty(app.lens.mLines{ind})
            lineName = strcat('Sampling Line','',num2str(ind));
            uitreenode(app.SamplingLinesNode,'Text',lineName,"NodeData",ind);
        end
    end
    expand(app.LineTree)
    app.ShowCropedIm.Value = true;
end
if ~isempty(app.lens.supLines)
    app.SupportLinesCheckBox.Value = true;
    app.SupportLinesCheckBox.Enable = true;
    app.DeleteLineButton.Enable = true;
    for ind = 1:length(app.lens.supLines)
        if ~isempty(app.lens.supLines{ind})
            lineName = strcat('Support Line','',num2str(ind));
            uitreenode(app.SupportLinesNode,'Text',lineName,"NodeData",ind);
        end
    end
    expand(app.LineTree)
end

if ~isempty(app.lens.calibIm)
    imshow(app.lens.calibIm,'Parent',app.UIAxesCalibration)
end
if ~isnan(app.lens.lambda)
    app.Wavelength.Value = app.lens.lambda;
end
if ~isnan(app.lens.lambda)
    app.n0.Value = app.lens.lambda;
end
if ~isnan(app.lens.calibKnownDistance)
    app.CalibDist.Value = app.lens.calibKnownDistance/10^-6;
end
if ~isempty(app.lens.calibLine)
    x = [app.lens.calibLine(1,1), app.lens.calibLine(2,1)];
    y = [app.lens.calibLine(1,2), app.lens.calibLine(2,2)];
    line(app.UIAxesCalibration,x,y,'linewidth',1.2,'color','b')
end
if ~isnan(app.lens.pixelSize)
    app.PixelSize.Value = app.lens.pixelSize;
end

end




