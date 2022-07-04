function selectPeaks(app)
app.addTwoConsole('Select peaks to follow. Press enter when done or "q" to quit')
app.listenForKey = true;
roi = drawrectangle(app.UIAxes);
uiwait(app.UIFigure)
if isequal(app.keyIn,'q')
    delete(roi)
    return
end
app.keyIn = '';
rec = roi.Position;
delete(roi);

lineInd = app.SamplingLineSpinner.Value;
layer = app.LayerSpinner.Value;
span = [rec(1) (rec(1)+rec(3))];
threshold = app.DistanceThresholdEditField.Value;
minPNums = app.MinPeakNumberEditField.Value;
flag = app.lens.mLines{lineInd}.AppTracePeaks(layer,span,threshold,minPNums);
if isequal(flag,'multiple hits')
    app.addTwoConsole('Error: The area created multiple hits, try smaller')
    return
elseif isequal(flag,'no hits')
    app.addTwoConsole('Error: There where no hits within the area')
    return
elseif isequal(flag,'allready in series')
    app.addTwoConsole('Error: The selected point already belongs to a serie')
    return
elseif isequal(flag,'not enough peaks')
    app.addTwoConsole('Error: Not enough peaks in the serie')
    return
end
app.addTwoConsole('New serie added')


% Recreate the table when a new series has been added.
% This is done in the AppTracePeaks function

% pinds = find(app.lens.mLines{lineInd}.Pinc);
% serie = flag;
% pointsInSerie = app.lens.mLines{lineInd}.S{serie}; % All the inds in that serie
% pN = length(pointsInSerie); % Number of points in serie
% group = ones(pN,1).*serie;
% pshift = zeros(pN,1);
% player = app.lens.mLines{lineInd}.PL(pointsInSerie);
% include = true(pN,1);
% 
% tab = {group, pointsInSerie, player, pshift, include};
% for tabind = 1:length(tab{1})
%     subtab = {tab{1}(tabind),tab{2}(tabind),tab{3}(tabind),...
%         tab{4}(tabind),tab{5}(tabind)};
%     app.PhaseTable.Data = [app.PhaseTable.Data;subtab];
% end

% app.estimatePhase() % Does not exist yet
end


% Old version

% function selectPeaks(app)
% app.addTwoConsole('Select peaks to follow. Press enter when done or "q" to quit')
% app.listenForKey = true;
% roi = drawrectangle(app.UIAxes);
% uiwait(app.UIFigure)
% if isequal(app.keyIn,'q')
%     delete(roi)
%     return
% end
% app.keyIn = '';
% rec = roi.Position;
% delete(roi);
%
% lineInd = app.SamplingLineSpinner.Value;
% layer = app.LayerSpinner.Value;
% p = app.lens.mLines{lineInd}.L{layer};
% pd = app.lens.mLines{lineInd}.PD(p);
% startP = rec(1) - app.lens.mLines{lineInd}.leftEdge;
% endP = startP + rec(3);
% key1 = pd > startP;
% key2 = pd < endP;
% key = key1 + key2;
% p = p(key == 2);
% if isempty(p)
%     app.addTwoConsole('Error: There wehere no peaks within the selected span')
%     return
% end
% pn = length(p);
% for pind = 1:pn
%     if app.lens.mLines{lineInd}.PS(p(pind)) == 0
%         continue
%     end
%     pSerie = app.lens.mLines{lineInd}.PS(p(pind)); % Serie the point belongs to
%     pointsInSerie = app.lens.mLines{lineInd}.S{pSerie}; % All the inds in that serie
%     pN = length(pointsInSerie); % Number of points in serie
%     group = ones(pN,1).*pSerie;
%     pshift = zeros(pN,1);
%     player = app.lens.mLines{lineInd}.PL(pointsInSerie);
%     interpType = ones(pN,1).*nan;
%     exclude = ones(pN,1);
%     exclude = exclude == 1;
%
%     tab = {group, pointsInSerie, pshift, player, interpType, exclude};
%     for tabind = 1:length(tab{1})
%         subtab = {tab{1}(tabind),tab{2}(tabind),tab{3}(tabind),...
%             tab{4}(tabind),tab{5}(tabind),tab{6}(tabind)};
%         app.PhaseTable.Data = [app.PhaseTable.Data;subtab];
%     end
% end
%
% end