function AppPlotLines(this,imhandle,lineType,lineId)
arguments
    this
    imhandle
    lineType
    lineId = nan;
end

switch lineType
    case 'support'
        if ~isempty(this.supLines)
            if ~isnan(lineId)
                x = [this.supLines{lineId}(1,1), this.supLines{lineId}(2,1)];
                y = [this.supLines{lineId}(1,2), this.supLines{lineId}(2,2)];
                line(imhandle,x,y,'linewidth',1.2,'color',[0.3010 0.7450 0.9330])
                text(imhandle,x(2), y(2), num2str(lineId),'FontSize',14,'Color','white');
            else
                for lineId = 1:length(this.supLines)
                    if ~isempty(this.supLines{lineId})
                        x = [this.supLines{lineId}(1,1), this.supLines{lineId}(2,1)];
                        y = [this.supLines{lineId}(1,2), this.supLines{lineId}(2,2)];
                        line(imhandle,x,y,'linewidth',1,'color','b')
                        text(imhandle,x(2), y(2), num2str(lineId),'FontSize',14,'Color','white');
                    end
                end
            end
        end
    case 'sampling'
        if ~isempty(this.mLines)
            if ~isnan(lineId)
                x = [this.mLines{lineId}.lineCoord(1,1), this.mLines{lineId}.lineCoord(2,1)];
                y = [this.mLines{lineId}.lineCoord(1,2), this.mLines{lineId}.lineCoord(2,2)];
                line(imhandle,x,y,'linewidth',1.2,'color',[0.9290 0.6940 0.1250])
                text(imhandle,x(2), y(2), num2str(lineId),'FontSize',14,'Color','white');
            else
                for lineId = 1:length(this.mLines)
                    if ~isempty(this.mLines{lineId})
                        x = [this.mLines{lineId}.lineCoord(1,1), this.mLines{lineId}.lineCoord(2,1)];
                        y = [this.mLines{lineId}.lineCoord(1,2), this.mLines{lineId}.lineCoord(2,2)];
                        line(imhandle,x,y,'linewidth',1,'color','r')
                        text(imhandle,x(2), y(2), num2str(lineId),'FontSize',14,'Color','white');
                    end
                end
            end
        end
end
end