% Function to calibrate size and other figure properties

function setFigure(span,varargin)
if ~isempty(varargin)
    if isequal(varargin{1},'subplot')
        hold on
        grid on
        ax = gca;
        ticks = getTickSize(span(1),span(2),15);
        ax.XTick = span(1):ticks:span(2);
    end
else
    set(gcf,'position',figureSize)
    hold on
    grid on
    ax = gca;
    ticks = getTickSize(span(1),span(2),15);
    ax.XTick = span(1):ticks:span(2);
end
end