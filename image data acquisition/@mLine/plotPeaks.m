%% Plot the extracted peaks
function plotPeaks(this, options)
arguments
    this    
    options.layerIndex uint8 = 1;
    options.interpolation {mustBeMember(options.interpolation,['raw','sgolay','gauss','all'])} = 'all';
end
if isempty(this.gaussPks)
    fprintf('Error: No peaks have been extracted \n')
end
layerIndex = options.layerIndex;

switch options.interpolation
    case 'sgolay'
        figure_title = strcat('Calculated peaks for layer : ',num2str(layerIndex),'',' (Savitzky–Golay interpolation)'); 
        figure('Name',figure_title,'NumberTitle','off')
        setFigure([0,this.d(end)])
        plot(this.d,this.sgolayPoints(:,layerIndex),'color',[0.6350 0.0780 0.1840])
        plot(this.sgolayPks{layerIndex}(:,2),this.sgolayPks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840].*0.8)
        this.plotSupportLines()
        
    case 'gauss'
        figure('Name','Calculated peaks (Gaussian interpolation)','NumberTitle','off')
        setFigure([0,this.d(end)])
        plot(this.d,this.gaussPoints(:,layerIndex),'color',[0 0.4470 0.7410].*0.8)
        plot(this.gaussPks{layerIndex}(:,2),this.gaussPks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.4)
        this.plotSupportLines()
        
    case 'all'
        figure('Name','Calculated peaks (All interpolations)','NumberTitle','off')
        setFigure([0,this.d(end)])
        plot(this.d,this.points(:,layerIndex),'*','color',[0.9290 0.6940 0.1250])
        plot(this.d,this.gaussPoints(:,layerIndex),'color',[0 0.4470 0.7410].*0.8)
        plot(this.d,this.sgolayPoints(:,layerIndex),'color',[0.6350 0.0780 0.1840])
%         plot(this.gaussPks{layerIndex}(:,2),this.gaussPks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.4)
%         plot(this.sgolayPks{layerIndex}(:,2),this.sgolayPks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840].*0.8)
        plot(this.pks{layerIndex}(:,2),this.pks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840].*0.8)
        this.plotSupportLines()
end
end