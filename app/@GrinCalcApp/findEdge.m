function findEdge(app)
lineId = app.SamplingLineSpinner.Value;

if ~app.ShowalllayersCheckBox.Value
    app.ShowalllayersCheckBox.Value = true;
end
app.uppdateLinePlot()

app.addTwoConsole('Select starting point of the span where to search for the lens edge')
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

gradPoints = zeros(size(app.lens.mLines{lineId}.points));
for layerInd = 1:app.lens.mLines{lineId}.imNums
    grad_temp = gradient(app.lens.mLines{lineId}.gaussPoints(:,layerInd)); % First gradient
    grad_temp = abs(gradient(grad_temp)); % Abs value of second gradient
    gradPoints(:,layerInd) = grad_temp;
end


leftEdgeList = zeros(app.lens.mLines{lineId}.imNums,1);
for layerInd = 1:app.lens.mLines{lineId}.imNums
    filtered_points = gradPoints(:,layerInd);
    filtered_points(app.lens.mLines{lineId}.d < xs) = 0;
    filtered_points(app.lens.mLines{lineId}.d > xe) = 0;
    [~, edgeInd] = max(filtered_points);
    leftEdgeList(layerInd) = edgeInd;
end
app.lens.mLines{lineId}.leftEdge = mean(leftEdgeList);
distToEdge = (app.lens.mLines{lineId}.d-app.lens.mLines{lineId}.leftEdge);
app.lens.mLines{lineId}.leftEdgeIndex = find(distToEdge > 0,1)-1;
app.signalFromApp = true;
app.LineDataTable.Data{2,2} = app.lens.mLines{lineId}.leftEdge;
% app.EdgeDistEditField.Value = app.lens.mLines{lineId}.leftEdge;
end