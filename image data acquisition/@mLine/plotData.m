function plotData(this,options)
arguments
    this
    options.interpolation {mustBeMember(options.interpolation,['raw','sgolay','gauss','all'])} = 'all'
    options.layerIndex = 'all';
end
if isempty(this.points)
    fprintf('ERROR: There are no sampled points \n')
    return
end

if isequal(options.interpolation,'all')
    if isequal(options.layerIndex, 'all')
        op = 1;
    elseif isnumeric(options.layerIndex)
        op = 2;
    end
elseif ~isequal(options.interpolation,'all')
    if isequal(options.layerIndex, 'all')
        op = 3;
    elseif isnumeric(options.layerIndex)
        op = 4;
    end
end
    
    
switch op
    case 1 % all_interpolations & all layers
        figure('Name','Plot sampled line (all layers)','NumberTitle','off')
        set(gcf,'position',[300,100,1000,800])
        subplot(3,1,1)
        setFigure([0,this.d(end)],'subplot')
        for ind = 1:this.imNums
            plot(this.d,this.points(:,ind),'color',[0.9290 0.6940 0.1250])
        end
        this.plotSupportLines()
        subplot(3,1,2)
        setFigure([0,this.d(end)],'subplot')
        for ind = 1:this.imNums
            plot(this.d,this.gaussPoints(:,ind),'color',[0 0.4470 0.7410].*0.8)
        end
        this.plotSupportLines()
        subplot(3,1,3)
        setFigure([0,this.d(end)],'subplot')
        for ind = 1:this.imNums
            plot(this.d,this.sgolayPoints(:,ind),'color',[0.6350 0.0780 0.1840])
        end
        this.plotSupportLines()
        
    case 2 % all_interpolations & single layer
        figure('Name',strcat('Plot sampled line (layer :',num2str(options.layerIndex),')'),'NumberTitle','off')
        set(gcf,'position',[300,100,1000,800])
        subplot(3,1,1)
        setFigure([0,this.d(end)],'subplot')
        plot(this.d,this.points(:,options.layerIndex),'color',[0.9290 0.6940 0.1250])
        this.plotSupportLines()
        subplot(3,1,2)
        setFigure([0,this.d(end)],'subplot')
        plot(this.d,this.gaussPoints(:,options.layerIndex),'color',[0 0.4470 0.7410].*0.8)
        this.plotSupportLines()
        subplot(3,1,3)
        setFigure([0,this.d(end)],'subplot')
        plot(this.d,this.sgolayPoints(:,options.layerIndex),'color',[0.6350 0.0780 0.1840])
        this.plotSupportLines()
        
    case 3 % Single interpolation & all layers
        figure_title = strcat('Plot sampled line (',options.interpolation,'',',all layers');
        figure('Name',figure_title,'NumberTitle','off')
        set(gcf,'position',[300,100,1000,800])
        setFigure([0,this.d(end)])
        this.plotSupportLines()
        switch options.interpolation
            case 'raw'
                for ind = 1:this.imNums
                    plot(this.d,this.points(:,ind),'color',[0.9290 0.6940 0.1250])
                end
            case 'gauss'
                for ind = 1:this.imNums
                    plot(this.d,this.gaussPoints(:,ind),'color',[0.9290 0.6940 0.1250])
                end
            case 'sgolay'
                for ind = 1:this.imNums
                    plot(this.d,this.sgolayPoints(:,ind),'color',[0.6350 0.0780 0.1840])
                end
        end
        
    case 4 % Single interpolation & single layer
        figure_title = strcat('Plot sampled line (',options.interpolation,'',',layer :',num2str(options.layerIndex),')');
        figure('Name',figure_title,'NumberTitle','off')
        set(gcf,'position',[300,100,1000,800])
        setFigure([0,this.d(end)])
        this.plotSupportLines()
        switch options.interpolation
            case 'raw'
                plot(this.d,this.points(:,options.layerIndex),'color',[0.9290 0.6940 0.1250])
            case 'gauss'
                plot(this.d,this.gaussPoints(:,options.layerIndex),'color',[0.9290 0.6940 0.1250])
            case 'sgolay'
                plot(this.d,this.sgolayPoints(:,options.layerIndex),'color',[0.6350 0.0780 0.1840])
        end
end
end