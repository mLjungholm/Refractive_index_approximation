% Delete lines from the stack class and the tables in the app
function deleteLine(app, style)
switch style
    case 'sampling'
        % Remove line from the Stack Class
        app.lens.AppDeleteLine('sampling',app.SamplingLineSpinner.Value)
        % Unenable the line spinner if ther is only one line left
        if app.lens.lineNums < 2
            app.SamplingLineSpinner.Enable = false;
        end
        % Find the table indices for the removed line in the table
        tableInd = find(cell2mat(app.SamplingLineTable.Data(:,1)) == app.SamplingLineSpinner.Value);
        % remove the line from the table
        app.SamplingLineTable.Data(tableInd,:) = [];
        % Set the line spinner to a new line if any left
        if isnan(app.lens.lineNums)            
            app.SamplingLineSpinner.Value = 1;
        else
            % Find the indice of the first of the remaining line/lines
            app.SamplingLineSpinner.Value = find(~cellfun(@isempty,app.lens.mLines),1);
        end        
        
    case 'support'
        % Remove line from Stack Class
        app.lens.AppDeleteLine('support',app.supLineInd)
        % Remove from table
        tableInd = find(app.SupportLineTable.Data{:,1} == app.supLineInd);
        app.upportLineTable.Data(tableInd,:) = [];
        % Find if there are any left and set the current line as the first
        % indice
        if ~isempty(app.upportLineTable.Data)
            app.supLineInd = app.upportLineTable.Data{1,1};
        else
            app.supLineInd = nan;
        end
end

% Uppdate the image
app.uppdateImstackImage()
end