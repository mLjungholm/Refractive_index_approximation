%% estimate the center of the line
function estimate_center(this)
this.plotData(interpolation = 'gauss')
fprintf('Set start of span /n')
[xs,~] = getXfromPlot;
fprintf('Set end of span /n')
[xe,~] = getXfromPlot;


gaussPeakList = zeros(this.imNums,1);
for layerInd = 1:this.imNums
    peakList = this.gaussPks{layerInd}(:,2);
    peakList(peakList < xs) = nan;
    peakList(peakList > xe) = nan;
    if nnz(~isnan(peakList)) ~= 1
        gaussPeakList(layerInd) = nan;
    else
        gaussPeakList(layerInd) = peakList(~isnan(peakList));
    end
end
gaussPeakList = gaussPeakList(~isnan(gaussPeakList));
this.centerLine = mean(gaussPeakList);
distToCenter = this.d-this.centerLine;
this.centerLineIndex = find(distToCenter > 0,1);
end
% function estimate_center(this,type)
% this.plotData(interpolation = type)
% fprintf('Set start of span /n')
% [xs,~] = getXfromPlot;
% fprintf('Set end of span /n')
% [xe,~] = getXfromPlot;
% if isequal(type,'gauss')
%     this.centerLine = getCenterofPoints(this.gaussPks(:,2),[xs,xe]);
% elseif isequal(type,'sgolay')
%     this.centerLine = getCenterofPoints(this.sgolayPks(:,2),[xs,xe]);
% end
% close
% end