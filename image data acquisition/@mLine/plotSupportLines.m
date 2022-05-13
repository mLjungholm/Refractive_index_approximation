%% Plotting support lines if any
function plotSupportLines(this)
if ~isnan(this.centerLine)
%     yspan = get(gcf,'CurrentAxes').YLim;
    plot([1;1].*this.centerLine,[0;1],'k')
end
if ~isnan(this.leftEdge)
%     yspan = get(gcf,'CurrentAxes').YLim;
    plot([1;1].*this.leftEdge,[0;1],'k')
end
if ~isnan(this.rightEdge)
%     yspan = get(gcf,'CurrentAxes').YLim;
    plot([1;1].*this.rightEdge,[0;1],'k')
end
end