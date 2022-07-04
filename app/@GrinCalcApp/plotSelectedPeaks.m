function plotSelectedPeaks(app)
colorlist = [0.3010 0.7450 0.9330;...
    0.4660 0.6740 0.1880;...
    0.4940 0.1840 0.5560;...
    0.8500 0.3250 0.0980;...
    0 0.4470 0.7410];
layerInd = app.LayerSpinner.Value;
lineInd = app.SamplingLineSpinner.Value;

if app.PeakSelectButton.Value
    key = app.lens.mLines{lineInd}.Pinc & app.lens.mLines{lineInd}.PL == layerInd; 

    % get the Distance values for the relevant peaks
    peaksD = app.lens.mLines{lineInd}.PD(key) + app.lens.mLines{lineInd}.leftEdge;
    % Get peak value for the relevant peaks
    peaksV = app.lens.mLines{lineInd}.PV(key);
    peaksS = app.lens.mLines{lineInd}.PS(key);
%     % Plot the peaks in the layer
    for pInd = 1:length(peaksD)
        plot(app.UIAxes,peaksD(pInd),peaksV(pInd),'color',colorlist(peaksS(pInd),:)...
            ,'Marker','s','MarkerSize',20,'LineWidth',5)
    end
    
elseif app.PeakTraceButton.Value
    key = app.lens.mLines{lineInd}.Pinc;
    peaksD = app.lens.mLines{lineInd}.PD(key) + app.lens.mLines{lineInd}.leftEdge;
    % Get peak value for the relevant peaks
    peaksV = app.lens.mLines{lineInd}.PV(key);
    peaksS = app.lens.mLines{lineInd}.PS(key);
    peaksL = app.lens.mLines{lineInd}.PL(key);
    for pInd = 1:length(peaksD)
        plot(app.UIAxes,peaksD(pInd),peaksV(pInd),'color',colorlist(peaksS(pInd),:)...
            ,'Marker','s','MarkerSize',20,'LineWidth',5)
        text(app.UIAxes,peaksD(pInd)+3,peaksV(pInd),num2str(peaksL(pInd)),...
            'HorizontalAlignment','left','FontSize',16)
    end
end
end

% if app.PeakSelectButton.Value
%     key1 = cell2mat(app.PhaseTable.Data(:,3)) == layerInd;
%     key2 = cell2mat(app.PhaseTable.Data(:,5)) == true;
%     key = key1 + key2;
%     key = key == 2;
%     % Sort out only peaks in the correct layer that are being included
%     peaks = cell2mat(app.PhaseTable.Data(key,1:2)); % [group, P_index];
%     if isempty(peaks)
%         return
%     end
%     % get the Distance values for the relevant peaks
%     peaksD = app.lens.mLines{app.SamplingLineSpinner.Value}.PD(peaks(:,2));
%     peaksD = peaksD + app.lens.mLines{app.SamplingLineSpinner.Value}.leftEdge;
%     % Get peak value for the relevant peaks
%     peaksV = app.lens.mLines{app.SamplingLineSpinner.Value}.PV(peaks(:,2));
%     peaksS = app.lens.mLines{app.SamplingLineSpinner.Value}.PS(peaks(:,2));
%     % Plot the peaks in the layer
%     for pInd = 1:size(peaks,1)
%         plot(app.UIAxes,peaksD(pInd),peaksV(pInd),'color',colorlist(peaksS(pInd),:)...
%             ,'Marker','s','MarkerSize',20,'LineWidth',5)
%     end
%     
% elseif app.PeakTraceButton.Value
%     key = app.PhaseTable.Data{:,5} == true;
%     peaks = app.PhaseTable.Data{key,1:2}; % [group, P_index];
%     % Get the uniqe groups in the layer
%     if isempty(peaks)
%         return
%     end
%     peaksD = app.lens.mLines{app.SamplingLineSpinner.Value}.PD(peaks(:,2));
%     peaksD = peaksD + app.lens.mLines{app.SamplingLineSpinner.Value}.leftEdge;
%     % Get peak value for the relevant peaks
%     peaksV = app.lens.mLines{app.SamplingLineSpinner.Value}.PV(peaks(:,2));
%     groups = unique(peaks(:,1));
%     groupsN = length(groups);
%     peaksL = app.lens.mLines{app.SamplingLineSpinner.Value}.PL(peaks(:,2));
%     for groupInd = 1:groupsN
%         gind = groups(groupInd);
%         pd = peaksD(peaks(:,1) == gind);
%         pv = peaksV(peaks(:,1) == gind);
%         scatter(app.UIAxes,pd,pv,'Marker','s','LineWidth',3,'SizeData',200,...
%             'CData',colorlist(gind,:),'MarkerFaceColor','none')
%     end
%     for pind = 1:size(peaks,1)
%         text(app.UIAxes,peaksD(pind)+3,peaksV(pind),num2str(peaksL(pind)),...
%             'HorizontalAlignment','left','FontSize',16)
%         
%     end



