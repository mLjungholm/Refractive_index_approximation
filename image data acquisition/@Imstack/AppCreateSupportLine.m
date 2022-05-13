% Creates a support line. (just used for guiding)
function ind = AppCreateSupportLine(this,roi)
% center_line = drawline(imhandle);
temp_coords = roi.Position;

if isempty(this.supLines)
    this.supLines = cell(1,1);
    this.supLines{1} = temp_coords; % Creates the line
    ind = 1;
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
end