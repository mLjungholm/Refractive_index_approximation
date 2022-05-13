function  getPhaseShift_Dev(this)

% tempLayerNum = 40; % Repladce with layerInd for the loop
for layerInd = 1:this.imNums
% for layerInd = 2:2


pointListRaw = this.points(:,layerInd);
pointList = this.sgolayPoints(:,layerInd);
d = this.d;
cropCenter;
pksList = getPeaks();
pksNum = uint8(size(pksList,1));
fprintf('Number of peaks : %u  Layer: %u \n ',pksNum, uint8(layerInd))


extraD = this.d(1:round(length(this.d)*0.6));
pointsRawExtra = this.points(1:length(extraD),layerInd);
plotPeaksDev()

roiLine = drawline;          % Calls draw function
pause                        % Waits for key press confirmation
tempLine = roiLine.Position; % Store the coordinates form the temporary roi object
close();
this.pks{layerInd} = getPksFromInterwal(tempLine(:,1),pksList(:,2));
pks = this.pks{layerInd};
% fprintf('Before double peak check. Peak nums: %u \n',uint8(length(pks)))
% disp('Peak coords:')
% disp(pks)
this.pks{layerInd} = checkDoubles(pks,30);
% fprintf('After double peak check. Peak nums: %u \n',uint8(length(pks)))
% disp('Peak coords:')
% disp(pks)


end

    function plotPeaksDev()
        if isempty(this.gaussPks)
            fprintf('Error: No peaks have been extracted \n')
        end
        figure('Name','Plot sampled lines','NumberTitle','off')
        setFigure([0,d(end)])
        plot(extraD,pointsRawExtra,'*','color',[0.9290 0.6940 0.1250])
        plot(d,pointList,'color',[0.6350 0.0780 0.1840])
        plot(pksList(:,2),pksList(:,1),'ob')
        this.plotSupportLines('croped')
    end



function pks = checkDoubles(pksList,threshold)
    pks = zeros(length(pksList),1).*nan;
%     pksI = 1;
    for pInd = 1:length(pksList)-1
        p1 = pksList(pInd);
        p2 = pksList(pInd+1);
        delta = abs(p2-p1);
        if delta < threshold
            pn = (p1+p2)/2;
            pksList(pInd+1) = pn;
        else
            pks(pInd) = p1;
%             pksI = pksI + 1;
        end        
    end
    pks(end) = pksList(end);
    pks = pks(~isnan(pks));
    
end

    
    function pks = getPksFromInterwal(span,pksList)
        map = pksList;
        map(map < span(1)) = 0;
        map(map > span(2)) = 0;
        pks = map(map ~= 0);
    end

    function cropCenter()
        d = d(d < this.centerLine);
        endInd = length(d);
        pointListRaw = pointListRaw(1:endInd,:);
        pointList = pointList(1:endInd,:);
    end

    function peaksS =  getPeaks()
            [pksH,locsH] = findpeaks(pointList,'MinPeakDistance',40);
            [~,locsL] = findpeaks(-pointList,'MinPeakDistance',40);
            peaksS = [pksH locsH; pointListRaw(locsL) locsL];
            peaksS = sortrows(peaksS,2);
            peaksS(:,2) = d(peaksS(:,2));
    end

end

