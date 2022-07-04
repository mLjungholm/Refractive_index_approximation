function estimatePhase(app,serie)
% Line index
lineInd = app.SamplingLineSpinner.Value;

% All points to be checked
pinds = app.lens.mLines{lineInd}.S{serie};

% Check start and end layer
leftMin = app.lens.mLines{lineInd}.leftPhaseMin;
leftMax = app.lens.mLines{lineInd}.leftPhaseMax;
if leftMin > leftMax
    startLayer = leftMax;
    endLayer = leftMin;
%     startLow = 0;
else
    startLayer = leftMin;
    endLayer = leftMax;
%     startLow = 1;
end

% Get staring phase from the defined peak (if there is a peak in the series
% in the starting layer)
foundInStart = false;
pind = ismember(app.lens.mLines{lineInd}.S{serie},app.lens.mLines{lineInd}.L{startLayer});
pN = nnz(pind);
if pN > 1
    app.addTwoConsole('Error: more than one point in the series was found in the starting layer')
    return
end
if pN
    pind = app.lens.mLines{lineInd}.S{serie}(pind);
    startPhase = app.lens.mLines{lineInd}.PP(pind);
    foundInStart = true;
end

if foundInStart && startPhase ~= 0
    lstart = double(startLayer);
    % Check each point in that serie
    for pi = 1:length(pinds)
        % If the point alreade has a phase value then ignore it.
        if app.lens.mLines{lineInd}.PP(pinds(pi)) == 0
            % Layer index of the point to be checked
            lind = double(app.lens.mLines{lineInd}.PL(pinds(pi)));
            % Estimated phase
            phase = startPhase + (lind-lstart)*app.LayerPhaseShift.Value;
            % Add value to the point
            app.lens.mLines{lineInd}.PP(pinds(pi)) = phase;
        end
    end
end


% Give phase to all points in end layer
%! This should be modified so that it can deal with the posibility that a
% peak series from the start layer also is in the end layer.
pind = ismember(app.lens.mLines{lineInd}.S{serie},app.lens.mLines{lineInd}.L{endLayer});
pN = nnz(pind);
foundInEnd = false;
if pN > 1
    app.addTwoConsole('Error: more than one point in the series was found in the end layer')
    return
end
if pN && foundInStart && startPhase ~= 0
    pind = app.lens.mLines{lineInd}.S{serie}(pind);
%     endPhase = app.lens.mLines{lineInd}.PP(pind);
    msg = strcat('Series belongs to both start and end layer. endPhase:'...
        ,num2str(app.lens.mLines{lineInd}.PP(pind)));
    app.addTwoConsole(msg)
    return
end
if pN
    pind = app.lens.mLines{lineInd}.S{serie}(pind);
    startPhase = app.lens.mLines{lineInd}.PP(pind);
    foundInEnd = true;
end

if foundInEnd
    lstart = double(endLayer);
    % Check each point in that serie
    for pi = 1:length(pinds)
        % If the point already has a phase value then ignore it.
        if app.lens.mLines{lineInd}.PP(pinds(pi)) == 0
            % Layer index of the point to be checked
            lind = double(app.lens.mLines{lineInd}.PL(pinds(pi)));
            % Estimated phase
            phase = startPhase + (lind-lstart)*app.LayerPhaseShift.Value;
            % Add value to the point
            app.lens.mLines{lineInd}.PP(pinds(pi)) = phase;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V2
% lineInd = SamplingLineSpinner.Value;
% 
% tabData = app.PhaseTable.Data;
% includeKey = cell2mat(tabData(:,5));
% tabData = cell2mat(tabData(:,1:4));
% 
% % All points to be checked
% pinds = tabData(includeKey,2);
% 
% leftMin = app.lens.mLines{lineInd}.leftPhaseMin;
% leftMax = app.lens.mLines{lineInd}.leftPhaseMax;
% if leftMin > leftMax
%     startLayer = leftMax;
%     endLayer = leftMin;
%     startLow = 0;
% else
%     startLayer = leftMin;
%     endLayer = leftMax;
%     startLow = 1;
% end
% 
% % Points in starting layer
% L = app.lens.mLines{lineInd}.L{startLayer};
% startKey = tabData(:,3) == startLayer;
% pindsStart = tabData(startKey,3);
% if isempty(L) || isempty(pindsStart)
%     app.addTwoConsole('Error: No peaks in the starting layer')
%     return
% end
% 
% % Give phase to all points in start layer
% for piStart = 1:length(pindsStart)
%     if app.lens.mLines{lineInd}.PS(pindsStart(piStart)) == 0
%         pOrder = find(L == pindsStart(piStart));
%         phase = app.StartHalfShift.Value + 1/2*(pOrder-1);
%         app.lens.mLines{lineInd}.PP(pindsStart(piStart)) = phase;
%     end
% end
% 
% % Give phase to all points in the same series as the start layer
% for piStart = 1:length(pindsStart)
%     % Get index of the serie
%     serieInd = app.lens.mLines{lineInd}.PS(pindsStart(piStart));
%     % All points from the table in that serie
%     pindsSerie = app.lens.mLines{lineInd}.S{serieInd};
%     % Layer of the starting point (shlud be the same as the startLayer
%     lstart = app.lens.mLines{lineInd}.PL(pindsStart(piStart));
%     % Phase of the starting point
%     phstart = app.lens.mLines{lineInd}.PP(pindsStart(piStart));
%     
%     % Check each point in that serie
%     for piSerie = 1:length(pindsSerie)
%         % If the point alreade has a phase value then ignore it.
%         if app.lens.mLines{lineInd}.PS(pindsSerie(piSerie)) == 0
%             % Layer index of the point to be checked
%             lind = app.lens.mLines{lineInd}.PL(pindsSerie(piSerie));
%             % Estimated phase
%             phase = phstart + (lind-lstart)*app.LayerPhaseShift.Value;
%             % Add value to the point
%             app.lens.mLines{lineInd}.PP(pindsSerie(piSerie)) = phase;
%         end
%     end
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%% V1
% function estimatePhase(app)
% tabData = app.PhaseTable.Data;
% includeKey = cell2mat(tabData(:,5));
% tabData = cell2mat(tabData(:,1:4));
% 
% leftMin = app.lens.mLines{SamplingLineSpinner.Value}.leftPhaseMin;
% leftMax = app.lens.mLines{SamplingLineSpinner.Value}.leftPhaseMax;
% if leftMin > leftMax
%     startLayer = leftMax;
%     endLayer = leftMin;
%     startLow = 0;
% else
%     startLayer = leftMin;
%     endLayer = leftMax;
%     startLow = 1;
% end
% 
% spanKey = tabData(:,3) >= starLayer & tabData(:,3) <= endLayer;
% if ~any(spanKey)
%     app.addTwoConsole('Error: No peaks where in the relevant layer span')
%     return
% end
% 
% % Find the points that are in hte staring layer
% startKey = tabData(:,3) == startLayer;
% if ~any(startKey)
%     app.addTwoConsole('Error: No peaks in the starting layer')
%     return
% end
% % Do a function that makes shure the first peak is 1/2 wavelength sift from
% % the edge. 
% key = includeKey & spanKey & startKey;
% if ~any(key)
%     app.addTwoConsole('Error: No peaks match criteria')
%     return
% end
% for pind = length(key)
%     if key(pind)
%         phase = 
%     end 
% end
