function AppPlotData(this,imhandle,interpolation,layerIndex)


switch interpolation
    case 'raw'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.d,this.points(:,ind),'color',[0.9290 0.6940 0.1250])
            end
        else
            plot(imhandle,this.d,this.points(:,layerIndex),'color',[0.9290 0.6940 0.1250])
        end
    case 'gauss'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.d,this.gaussPoints(:,ind),'color',[0.9290 0.6940 0.1250])
            end
        else
            plot(imhandle,this.d,this.gaussPoints(:,layerIndex),'color',[0.9290 0.6940 0.1250])
        end
    case 'sgolay'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.d,this.sgolayPoints(:,ind),'color',[0.6350 0.0780 0.1840])
            end
        else
            plot(imhandle,this.d,this.sgolayPoints(:,layerIndex),'color',[0.6350 0.0780 0.1840])
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
else
    ymax = max(max(this.sgolayPoints));
end

yt = linspace(0,ymax,10);
imhandle.YTick = yt;
end




