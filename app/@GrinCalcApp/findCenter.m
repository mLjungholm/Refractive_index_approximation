function findCenter(app)
lineId = app.SamplingLineSpinner.Value;

if ~app.ShowalllayersCheckBox.Value
    app.ShowalllayersCheckBox.Value = true;
end
app.uppdateLinePlot()



app.addTwoConsole('Select starting point of the span where to search for the lens center')
app.addTwoConsole('Move the cursor to the starting point. Press "Enter" when done')
app.listenForKey = true;
h1 = drawcrosshair('Parent',app.UIAxes,'Position',[50 50],'LineWidth',1,'Color','k');
uiwait(app.UIFigure)
xs = h1.Position(1);

if isequal(app.keyIn,'q')
    delete(h1)
    return
end
delete(h1)
yt = app.UIAxes.YTick(end);
plot(app.UIAxes,[xs;xs],[0;yt],'k')

app.keyIn = '';
app.addTwoConsole('Move the cursor to the End point. Press "Enter" when done')
app.listenForKey = true;
h2 = drawcrosshair('Parent',app.UIAxes,'Position',[50 50],'LineWidth',1,'Color','k');
uiwait(app.UIFigure)
xe = h2.Position(1);
if isequal(app.keyIn,'q')
    delete(h2)
    return
end
delete(h2)
app.keyIn = '';
plot(app.UIAxes,[xe;xe],[0;yt],'k')

[~,is] = min(abs(app.lens.mLines{lineId}.d - xs));
[~,ie] = min(abs(app.lens.mLines{lineId}.d - xe));

centerPoints = zeros(app.lens.mLines{lineId}.imNums,1);
for layerInd = 1:app.lens.mLines{lineId}.imNums
    centerSpan = app.lens.mLines{lineId}.gaussPoints(is:ie,layerInd);
    spanGrad = gradient(centerSpan);
    if nnz(spanGrad >= 0) == length(centerSpan) || nnz(spanGrad <= 0) == length(centerSpan)
        centerPoints(layerInd) = nan;
        continue
    end
    spanGrad = abs(spanGrad);    
    [~, centerInd] = min(spanGrad);
    centerPoints(layerInd) = centerInd + xs-1;

end
centerPoints = centerPoints(~isnan(centerPoints));
app.lens.mLines{lineId}.centerLine = mean(centerPoints);
distToCenter = app.lens.mLines{lineId}.d-app.lens.mLines{lineId}.centerLine;
app.lens.mLines{lineId}.centerLineIndex = find(distToCenter > 0,1);

app.signalFromApp = true;
app.LineDataTable.Data{1,2} = app.lens.mLines{lineId}.centerLine;
% app.CenterDistEditField.Value = app.lens.mLines{lineId}.centerLine;
end