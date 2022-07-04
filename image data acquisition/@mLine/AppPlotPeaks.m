function AppPlotPeaks(this,imhandle,interpolation,layerIndex)
switch interpolation
    case 'raw'
        return
    case 'gauss'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.gaussPks{ind}(:,2),this.gaussPks{ind}(:,1),'o','color',[0 0.4470 0.7410].*0.8)
            end
        else
            plot(imhandle,this.gaussPks{layerIndex}(:,2),this.gaussPks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.8)
        end
    case 'sgolay'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.sgolayPks{ind}(:,2),this.sgolayPks{ind}(:,1),'o','color',[0.6350 0.0780 0.1840])
            end
        else
            plot(imhandle,this.sgolayPks{layerIndex}(:,2),this.sgolayPks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840])
        end
    case 'all'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.gaussPks{ind}(:,2),this.gaussPks{ind}(:,1),'o','color',[0 0.4470 0.7410].*0.8)
                plot(imhandle,this.sgolayPks{ind}(:,2),this.sgolayPks{ind}(:,1),'o','color',[0.6350 0.0780 0.1840])
            end
        else
            plot(imhandle,this.sgolayPks{layerIndex}(:,2),this.sgolayPks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840])
            plot(imhandle,this.gaussPks{layerIndex}(:,2),this.gaussPks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.8)
        end
    case 'pks'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.pks{ind}(:,2),this.pks{ind}(:,1),'o','color',[0 0.4470 0.7410].*0.8)
            end
        else
            plot(imhandle,this.pks{layerIndex}(:,2),this.pks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.8)
        end
        
end
xmax = this.d(end);
xt = linspace(0,xmax,10);
xt = round(xt);
imhandle.XTick = xt;


if isequal(interpolation,'raw')
    ymax = max(max(this.points));
elseif isequal(interpolation,'gauss')
    ymax = max(max(this.gaussPoints));
elseif isequal(interpolation,'sgolay')
    ymax = max(max(this.sgolayPoints));
else
    ymax1 = max(max(this.sgolayPoints));
    ymax2 = max(max(this.gaussPoints));
    ymax3 = max(max(this.points));
    ymax = max([ymax1 ymax2 ymax3]);
end

yt = linspace(0,ymax,10);
imhandle.YTick = yt;

if ~isnan(this.leftEdge)
    plot(imhandle,[this.leftEdge;this.leftEdge],[0;yt(end)],'k')
end
if ~isnan(this.centerLine)
    plot(imhandle,[this.centerLine;this.centerLine],[0;yt(end)],'k')
end

end




