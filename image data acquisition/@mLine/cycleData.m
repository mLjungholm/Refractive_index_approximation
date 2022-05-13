% Cycle through the peak data
function cycleData(this)

figure('Name','Calculated peaks (All interpolations)','NumberTitle','off')
set(gcf,'position',[300,100,1000,800])

cspan = getCenterSpan;

for layerIndex = 1:this.imNums
    plot(this.d,this.gaussPoints(:,layerIndex),'color',[0 0.4470 0.7410].*0.8)
    hold on
    plot(this.d,this.points(:,layerIndex),'*','color',[0.9290 0.6940 0.1250])
    %  plot(this.gaussPks{layerIndex}(:,2),this.gaussPks{layerIndex}(:,1),'o','color',[0 0.4470 0.7410].*0.4)
    plot(this.d,this.sgolayPoints(:,layerIndex),'color',[0.6350 0.0780 0.1840])
    %  plot(this.sgolayPks{layerIndex}(:,2),this.sgolayPks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840].*0.8)
%     plot(this.pks{layerIndex}(:,2),this.pks{layerIndex}(:,1),'o','color',[0.6350 0.0780 0.1840].*0.8)
    lPeaks = this.L{layerIndex};
    plot(this.PD(lPeaks)+this.leftEdge,this.PV(lPeaks),'o','color',[0.6350 0.0780 0.1840].*0.8)
    
    txt = {'Layer Nr:',num2str(layerIndex)};
    text(200,0.5,txt)
    
    if ~isnan(this.sgolayZone)
        v = [this.leftEdge 0; this.leftEdge 1; this.sgolayZone 1; this.sgolayZone 0];
        f = [1 2 3 4];
        patch('Faces',f,'Vertices',v,'FaceColor','[0.8500 0.3250 0.0980]','FaceAlpha',0.2);
    end
    if ~isnan(this.centerZone)
        le = this.centerLine - this.centerZone;
        re = this.centerLine + this.centerZone;
        v = [le 0; le 1; re 1; re 0];
        f = [1 2 3 4];
        patch('Faces',f,'Vertices',v,'FaceColor','[0.8500 0.3250 0.0980]','FaceAlpha',0.2);
        text(this.centerLine-20,-0.025,'Center exclution zone')
        
    end
    endVal = cspan(layerIndex);
    plot([le;re],[endVal;endVal],'k')
    plot([1;1].*this.centerLine,[0;1],'k','linewidth',2)
    text(this.centerLine-20,1.025,'Estimated Center')
    plot([1;1].*this.leftEdge,[0;1],'k','linewidth',2)
    text(this.leftEdge-20,1.025,'Left Edge')
    plot([1;1].*this.rightEdge,[0;1],'k','linewidth',2)
    text(this.rightEdge-20,1.025,'Right Edge')
    
    ax = gca;
    ax.YLim = [-0.05 1.05];
    grid minor
    
    pause
    hold off
end
close

    function centerSpan = getCenterSpan()
        dwidth = 5;
        cspan = zeros(this.imNums,1);
        for layerInd = 1:this.imNums
            s = this.points(this.centerLineIndex-dwidth:this.centerLineIndex+dwidth,layerInd);
            s = mean(s);
            cspan(layerInd) = s;
        end
        cspan = cspan - min(cspan);
        centerSpan = cspan./max(cspan);
        this.centerVal = cspan;
    end
end