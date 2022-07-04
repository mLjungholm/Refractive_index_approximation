% Delete line from the stack
function AppDeleteLine(this,linetype,lineInd)
switch linetype
    case 'sampling'
        if this.lineNums > 1            
            this.mLines{lineInd} = []; % The removed line creates an empty spot.
            this.lineNums = this.lineNums - 1;
        else
            % Return mLines to empty list of the last line is removed
            this.lineNums = nan;
            this.mLines = [];
        end
        
    case 'support'
        this.supLines{lineInd} = [];
        emptyCell = cellfun(@isempty,this.supLines);
        if ~any(emptyCell)
            this.supLines = [];
        end
end
end