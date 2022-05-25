%% Function for extracting the peaks & valeys form the smoothed data (needs to be fail profed for empty data)
function findPeaks(this)
this.gaussPks = {};
this.sgolayPks = {};
if isnan(this.leftEdge)
    this.leftEdge = 0;
end
for layerInd = 1:this.imNums
    % Find the peaks and valeys of the gaussian smoothed curve
    [pksH,locsH] = findpeaks(this.gaussPoints(:,layerInd),'MinPeakDistance',40);
    [~,locsL] = findpeaks(-this.gaussPoints(:,layerInd),'MinPeakDistance',40);
    peaksG = [pksH locsH ; this.gaussPoints(locsL,layerInd) locsL ];
    peaksG = sortrows(peaksG,2);
    d = this.d(peaksG(:,2));
    cutmap = d > this.leftEdge;
    peaksG = peaksG(cutmap,:);
    peaksG(:,2) = this.d(peaksG(:,2));
    peaksG = checkDoubles(peaksG,30);
    this.gaussPks{layerInd} = peaksG;
    
    % Find the peaks and valeys of the Savitzky–Golay curve. (The gaussian
    % is better att larger intervalls wheras the Savitzky–Golay is better
    % close to the edges
    [pksH,locsH] = findpeaks(this.sgolayPoints(:,layerInd),'MinPeakDistance',40);
    [~,locsL] = findpeaks(-this.sgolayPoints(:,layerInd),'MinPeakDistance',40);
    peaksS = [pksH locsH; this.sgolayPoints(locsL,layerInd) locsL];
    peaksS = sortrows(peaksS,2);
    d = this.d(peaksS(:,2));
    cutmap = d > this.leftEdge;
    peaksS = peaksS(cutmap,:);
    peaksS(:,2) = this.d(peaksS(:,2));
    peaksS = checkDoubles(peaksS,30);
    this.sgolayPks{layerInd} = peaksS;
    
    this.pks{layerInd} = summarizePeaks(peaksG,peaksS);
end



function pks = checkDoubles(pksList,threshold)
    pks = zeros(size(pksList));
    for pInd = 1:size(pksList,1)-1
        p1 = pksList(pInd,2);
        p2 = pksList(pInd+1,2);
        delta = abs(p2-p1);
        if delta < threshold
            pn = (p1+p2)/2;
            pksList(pInd+1,2) = pn;
            pksList(pInd+1,1) = (pksList(pInd+1,1)+pksList(pInd,1))/2;
        else
            pks(pInd,2) = p1;
            pks(pInd,1) = pksList(pInd,1);
%             pksI = pksI + 1;
        end        
    end
    pks(end,:) = pksList(end,:);
    map = pks(:,2) ~= 0;
    pks = pks(map,:);
end


    function pks = summarizePeaks(gaussList,sgolayList)
        smap = sgolayList(:,2) < 150;
        spart = sgolayList(smap,:);
        gmap = gaussList(:,2) > 150;
        gpart = gaussList(gmap,:);
        pks = [spart;gpart];
    end
end
