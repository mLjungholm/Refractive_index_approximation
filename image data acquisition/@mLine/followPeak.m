%% In development (uncompleated)

% The idea was to make a function where you select a single peak that is to
% be traced througout the layers.

% This was abanoned for a more compleate function that tests all extracted
% peaks, function "mLine.tracePeaks()"

function followPeak(this,type,startLayer)
arguments
    this    
    type string {mustBeMember(type,['sgolay','gauss'])} = 'gauss';
    startLayer uint8 = 1
end

this.plotPeaks(layerIndex = startLayer, interpolation = type)
fprintf('Draw a rectangle around the peak that is to be traced. Press any key when done \n')
rec =  drawrectangle;
pause
span = [rec.Position(1) rec.Position(1)+rec.Position(3)];
stepSize = this.centerLine/this.imNums;

peakPos = zeros(this.imNums,1);
[flag, currentP] = findPeak(span);
if ~flag
    fprintf('Error: No peaks within the span \n')
    return
end
peakPos(1) = currentP;
for layerInd = startLayer:this.imNums
    [flag, currentP] = findPeak(span);
    if
end



    function [flag, p] = findPeak(span)
        
    end

end