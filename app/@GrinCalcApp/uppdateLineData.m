function uppdateLineData(app)
% if app.SamplingLineSpinner.Value == app.currentLineInd
%     return
% elseif app.SamplingLineSpinner.Value > length(app.lens.mLines)
%     app.SamplingLineSpinner.Value = app.currentLineInd;
%     return
% elseif isempty(app.lens.mLines{app.SamplingLineSpinner.Value})
%     app.SamplingLineSpinner.Value = app.currentLineInd;
%     return
% end

app.uppdateLinePlot()
app.currentLineInd = app.SamplingLineSpinner.Value;
% app.SamplingLineSpinnerPhase.Value = app.SamplingLineSpinner.Value;

if ~isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.leftEdge)
%     app.signalFromApp = true;
    app.LineDataTable.Data{2,2} = app.lens.mLines{app.SamplingLineSpinner.Value}.leftEdge;
%     app.EdgeDistEditField.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.leftEdge;
end


if ~isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.centerLine)
    app.signalFromApp = true;
    app.LineDataTable.Data{1,2} = app.lens.mLines{app.SamplingLineSpinner.Value}.centerLine;
end

if ~isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMin)
%     app.signalFromApp = true;
    app.LineDataTable.Data{3,2} = app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMin;
%     app.signalFromApp = true;
    app.LineDataTable.Data{4,2} = app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMax;
end

% if ~isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.lambda)
% app.Wavelengt_2.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.lambda;
% end
% if ~isnan(app.lens.mLines{app.SamplingLineSpinner.Value}.n0)
%     app.n0_2.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.n0;
% end

% if isempty(app.lens.mLines{app.SamplingLineSpinner.Value}.gaussPks) || app.ShowPeaks.Value
%     app.ShowPeaks.Value = false;
% elseif ~isempty(app.lens.mLines{app.SamplingLineSpinner.Value}.gaussPks) || ~app.ShowPeaks.Value
%     app.ShowPeaks.Value = true;
% end
    
end