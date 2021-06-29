% Overcomplicated (but it works) function for extracting the peaks and valeys of the phase
% shift values.

function [peakVal,peakPos,nrPeaks] = findPeaks(coords,vals,threshold)
inN = length(vals);

if coords(1) < coords(end)
    coords = flipud(coords);
    vals = flipud(vals);
end
peakVal = zeros(inN,1).*nan;
peakPos = zeros(inN,1).*nan;
peakVal(1) = vals(1); peakVal(end) = vals(end);
peakPos(1) = coords(1); peakPos(end) = coords(end);

yMax = max(vals);
yMin = min(vals);
yRange = yMax - yMin;

endoflist = 0;
ind0 = 2;
newListInd = 2;
while ~endoflist
    tempList = zeros(inN,2).*nan;
    v0 = vals(ind0);
    v1 = vals(ind0 +1);
    if v1 > v0 && v1 > (yMax - threshold*yRange) && v0 > (yMax - threshold*yRange)
        flag = 1;
    elseif v1 < v0 && v1 < (yMin + threshold*yRange) && v0 < (yMin + threshold*yRange)
        flag = -1;
    else
        flag = 0;
    end
    switch flag
        case 1
            tempList(ind0,:) = [v0 coords(ind0)];
            for ind = ind0+1:inN
                v1 = vals(ind);
                p0 = coords(ind);
                if v1 > (yMax - threshold*yRange)
                    tempList(ind,:) = [v1 p0];
                else
                    ind0 = ind;
                    break
                end
                if ind == inN
                    endoflist = 1;
                end
            end
            [v,pi] = max(tempList(:,1));
            peakVal(newListInd) = v;
            peakPos(newListInd) = tempList(pi,2);
            newListInd = newListInd + 1;
        case -1
            tempList(ind0,:) = [v0 coords(ind0)];
            for ind = ind0+1:inN
                v1 = vals(ind);
                p0 = coords(ind);
                if v1 < (yMin + threshold*yRange)
                    tempList(ind,:) = [v1 p0];
                else
                    ind0 = ind;
                    break
                end
                if ind == inN
                    endoflist = 1;
                end
            end
            [v,pi] = min(tempList(:,1));
            peakVal(newListInd) = v;
            peakPos(newListInd) = tempList(pi,2);
            newListInd = newListInd + 1;
        case 0
            ind0 = ind0 + 1;
            if ind0 == inN
                endoflist = 1;
            end
    end
end
peakVal = peakVal(~isnan(peakVal));
peakPos = peakPos(~isnan(peakPos));
nrPeaks = length(peakVal);
end