function AppDeleteLine(this,linetype,lineInd)
switch linetype
    case 'sampling'
        if this.lineNums > 1
            this.mLines{lineInd} = [];
            this.lineNums = this.lineNums - 1;
        else
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