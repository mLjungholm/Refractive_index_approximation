% Creates a sampling line and calls the functions for sampling the
% data and smoothing
function drawLine(this)

if isempty(this.imSize) % The function only works on a cropped image set
    fprintf('ERROR: The image stack has not been croped\n')
    return
end
fprintf('Draw line in the figure. Press any key when done to continue. \n')
figure('Name','Draw measurement line','NumberTitle','off')
imshow(this.imStack(:,:,1))
set(gcf,'position',this.figSize)
axis equal
if ~isempty(this.supLines) % Checks if there is a center line defined and draws it
    for ind = 1:length(this.supLines)
        if ~isempty(this.supLines{ind})
            x = [this.supLines{ind}(1,1), this.supLines{ind}(2,1)];
            y = [this.supLines{ind}(1,2), this.supLines{ind}(2,2)];
        end
        line(x,y,'linewidth',1,'color','b')
    end
end
if ~isnan(this.lineNums)
    for lineInd = 1:size(this.mLines,1)
        if isempty(this.mLines{lineInd})
            continue
        end
        this.mLines{lineInd}.plotLine
    end
end
roiLine = drawline;          % Calls draw function
pause                        % Waits for key press confirmation
tempLine = roiLine.Position; % Store the coordinates form the temporary roi object
close();

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