% Simply extracts the x,y coordinates from a plot
function [x,y] = getXfromPlot
h = drawcrosshair('Position',[50 50],'LineWidth',1,'Color','k');
fprintf('Drag the cross to the desired point. Press any key when done \n')
pause
x = h.Position(1);
y = h.Position(2);

end

% imStack.plotPeaks(1,1)
% h = drawcrosshair('Position',[50 50],'LineWidth',1,'Color','k');
% fprintf('Drag the cross to the desired point. Press any key when done \n')
% pause
% x = h.Position(1);
% axval = get(gcf,'CurrentAxes');
% ay = axval.YLim;

% close

% imStack.plotPeaks(1,1)
% plot([x;x],ay','k')
