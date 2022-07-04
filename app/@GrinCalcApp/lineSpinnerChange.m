function lineSpinnerChange(app)
% Check if the line id the spinner changed to exists
inds = cell2mat(app.SamplingLineTable.Data(:,1));

%D Does the new lower value exist?
if ~isempty(find(inds == app.SamplingLineSpinner.Value,1))
    app.currentLineInd = app.SamplingLineSpinner.Value;
    return
end

% is the new value higher or lower than the old valure?
if app.SamplingLineSpinner.Value > app.currentLineInd
    % Does any higher values exist?
    higherValInd = find(inds > app.SamplingLineSpinner.Value,1);
    if ~isempty(higherValInd)
            app.SamplingLineSpinner.Value = higherValInd;
            app.currentLineInd = higherValInd;
    else
            % No higher values --> return the sppinner to the old value
            app.SamplingLineSpinner.Value = app.currentLineInd;
    end
else
    % Do any lower values exist?
    lowerValueInd = find(inds < app.SamplingLineSpinner.Value,1,'last');
    if ~isempty(lowerValueInd)
            app.SamplingLineSpinner.Value = lowerValueInd;
            app.currentLineInd = lowerValueInd;
    else
            % No lower values --> return the sppinner to the old value
            app.SamplingLineSpinner.Value = app.currentLineInd;
    end
end
end