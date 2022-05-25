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
                plot(imhandle,this.d,this.gaussPoints(:,ind),'color',[0 0.4470 0.7410].*0.8)
            end
        else
            plot(imhandle,this.d,this.gaussPoints(:,layerIndex),'color',[0 0.4470 0.7410].*0.8)
        end
    case 'sgolay'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.d,this.sgolayPoints(:,ind),'color',[0.6350 0.0780 0.1840])
            end
        else
            plot(imhandle,this.d,this.sgolayPoints(:,layerIndex),'color',[0.6350 0.0780 0.1840])
        end
    case 'all'
        if isequal(layerIndex,'all')
            for ind = 1:this.imNums
                plot(imhandle,this.d,this.points(:,ind),'color',[0.9290 0.6940 0.1250])
                plot(imhandle,this.d,this.gaussPoints(:,ind),'color',[0 0.4470 0.7410].*0.8)
                plot(imhandle,this.d,this.sgolayPoints(:,ind),'color',[0.6350 0.0780 0.1840])
            end
        else
            plot(imhandle,this.d,this.points(:,layerIndex),'color',[0.9290 0.6940 0.1250])
            plot(imhandle,this.d,this.gaussPoints(:,layerIndex),'color',[0 0.4470 0.7410].*0.8)
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


if ~isnan(this.centerZone) && ~isnan(this.centerLine)
%     le = this.centerLine - this.centerZone;
%     re = this.centerLine + this.centerZone;
    v = [this.centerZone 0; this.centerZone ymax; this.centerLine ymax; this.centerLine 0];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor','[0.8500 0.3250 0.0980]','FaceAlpha',0.2,'Parent',imhandle);
%     text(this.centerLine-20,-0.025,'Center exclution zone')
    
end

end




