% Creates a support line. (just used for guiding)
function createSupportLine(this)
if isempty(this.imSize)
    fprintf('ERROR: The image stack has not been croped\n')
    return
end
fprintf('Draw centerLine in the figure. Press any key when done to continue. \n')
figure('Name','Draw center line','NumberTitle','off')
imshow(this.imStack(:,:,1))
set(gcf,'position',this.figSize)
axis equal
if ~isempty(this.supLines)
    for ind = 1:length(this.supLines)
        if ~isempty(this.supLines{ind})
            x = [this.supLines{ind}(1,1), this.supLines{ind}(2,1)];
            y = [this.supLines{ind}(1,2), this.supLines{ind}(2,2)];
            line(x,y,'linewidth',1,'color','b')
        end
    end
end

center_line = drawline;
temp_coords = center_line.Position;
pause

if isempty(this.supLines)
    this.supLines = cell(1,1);
    this.supLines{1} = temp_coords; % Creates the line
else
    emptyCell = cellfun(@isempty,this.supLines); % Search for any empty cells in mLines array
    if any(emptyCell)
        ind = find(emptyCell,1);
        this.supLines{ind,1} = temp_coords; % Place the new line in the empty cell
    else
        ind = length(this.supLines) + 1;
        this.supLines{ind,1} = temp_coords; % Create a new cell an store line
    end
end
close();
end