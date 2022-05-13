function id = AppCreateSamplingLine(this,roi)
tempLine = roi.Position; % Store the coordinates form the temporary roi object

% Check if there are any previous lines. Else creates the mLine
% cell array
if isnan(this.lineNums)
    this.mLines = cell(1,1);
    this.mLines{1} = mLine(tempLine,1); % Creates the line
    this.lineNums = 1;
    id = 1;
else
    emptyCell = cellfun(@isempty,this.mLines); % Search for any empty cells in mLines array
    if any(emptyCell)
        id = find(emptyCell);
        this.mLines{id,1} = mLine(tempLine,id); % Place the new line in the empty cell
        this.lineNums = this.lineNums + 1;
    else
        this.lineNums = this.lineNums + 1;
        id = this.lineNums;
        this.mLines{this.lineNums,1} = mLine(tempLine,this.lineNums); % Create a new cell an store line
    end
end
this.sampleData(id); % Call sampling function
this.lastLineId = id;
end