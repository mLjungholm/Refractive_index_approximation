function flag = AppTracePeaks(this, layer, span, threshold, min_nums_in_serie)
% If the peaks have not been clasified then do that first
if isempty(this.PD)
    pks = this.pks; % input peaks
    pks = cellfun(@(x) x(:,2),pks,'UniformOutput',false); % extract only distance
    for layerInd = 1:this.imNums
        tempPks = pks{layerInd};
        tempPks = tempPks(tempPks > this.leftEdge);
        pks{layerInd} = tempPks(tempPks < this.centerLine);
    end
    
    numsP = sum(cell2mat(cellfun(@length ,pks,'UniformOutput',false))); % Total number of peaks
    this.pksNums = numsP;
    numsL = this.imNums; % number of layers
%     numsS = 0; % number of series
    
    % P = cell(1,numsP); % Each cell in P(points) = [dist, layer, series, phase];
    % P = cellfun(@(x) zeros(4,1),P,'UniformOutput',false);
    this.PD = zeros(numsP,1);
    this.PL = uint16(zeros(numsP,1));
    this.PS = uint16(zeros(numsP,1));
    this.PP = zeros(numsP,1);
    this.PV = zeros(numsP,1); % peak value
    this.L = cell(1,numsL); % Each cell in this.L(layer) = [points];
    this.S = {}; % Repository for all peak series this.S(series) = [points];
    this.Pinc = false(numsP,1);
%     checkedPeaks = zeros(numsP,1);
    
    
    % Give indexes and dist values to "P" and "this.L" for all peaks
    startP = 1;
    for indL = 1:numsL
        tempList = zeros(length(pks{indL}),1);
        for i = 1:length(pks{indL}) % Loop through points in layer
            this.PD(startP)= pks{indL}(i);
            this.PV(startP) = this.pks{indL}(i,1);
            this.PL(startP) = indL;
            %         P{indP}(1) = pks{indL}(i);
            %         P{indP}(2) = indL;
            tempList(i) = startP;
            startP = startP + 1;
        end
        this.L{indL} = tempList;
    end
    % Correct the distances to the left edge
    this.PD = this.PD - this.leftEdge;
end
%%%%%%%%%%%%%%%%%%%%%%%
% Modify so that the first(above) part is only done once


% Find the peak (if any) within the input-span
peaks = this.L{layer}; % Peaks in the input layer
peaksD = this.PD(peaks); % dist values for the peaks
% starP = span(1); % start and end values for the span
% endP = span(2);

key1 = peaksD > (span(1)-this.leftEdge);
key2 = peaksD < (span(2)-this.leftEdge);
key = key1 + key2;
peakInd = peaks(key == 2); % This should result in a single value
if length(peakInd) > 1
    flag = 'multiple hits';
    return
elseif isempty(peakInd)
    flag = 'no hits';
    return
end

% Check if the peak belongs to any other series
searchSeries = cellfun(@(x) find(x==peakInd),this.S,'UniformOutput',false);
foundPeaks = ~cellfun(@isempty,searchSeries);
if any(foundPeaks)
    flag = 'allready in series';
    return
end

startP = peakInd;
numsL = this.imNums;
traceList = zeros(numsL,1); % Create empty list (one possible slot for each layer)
nums_in_serie = 1;  % Keep track of the numer of points in the series
start_layer = layer;   % What layer does the trace start at


% If the image sequence aw perfect then you should be able to loop
% around the last image and contiue on the first. In this case it is
% not so. So the code will end att the last or first layer.

% Trace forward layers
d0 = this.PD(startP); % Set the first point distance
traceList(layer) = peakInd;
for layer_ind = start_layer:numsL-1
    test_layer = layer_ind + 1; % Test layer
    points_in_layer = this.L{test_layer};    % all point indices in that layer
    dist_in_layer = this.PD(points_in_layer); % Corresponding distance for each point
    [deltaP, closest_ind] = min(abs(d0 - dist_in_layer));   % Find the point that is closest
    closest_point = points_in_layer(closest_ind); % Get the global index for that point
    if deltaP < threshold  % Check if the distance is below the threshold.
        traceList(test_layer) = closest_point; % Add point to list
        nums_in_serie = nums_in_serie + 1; % Increase count
        d0 = this.PD(closest_point);
    else
        break
    end
end

% Trace backwards through the layers from the starting point
d0 = this.PD(startP); % Set the first point distance
for layer_ind = start_layer:-1:2
    test_layer = layer_ind - 1; % Test layer
    points_in_layer = this.L{test_layer};    % all point indices in that layer
    dist_in_layer = this.PD(points_in_layer);    % Corresponding distance for each point
    [deltaP, closest_ind] = min(abs(d0 - dist_in_layer));   % Find the point that is closest
    closest_point = points_in_layer(closest_ind); % Get the global index for that point
    if deltaP < threshold % Check if the distance is below the threshold
        traceList(test_layer) = closest_point; % Add point to list
        nums_in_serie = nums_in_serie + 1; % Increase count
        d0 = this.PD(closest_point);
    else
        break
    end
end

% Create the new serie
if nums_in_serie >= min_nums_in_serie
    if isempty(this.S)
        numsS = 1;
    else
        numsS = length(this.S) + 1;        
    end
    this.S{numsS} = traceList(traceList ~= 0);
else
    flag = 'not enough peaks';
    return
end
flag = numsS;
% Add what series each point belongs to
pinds = this.S{numsS};
for pi = 1:length(pinds)
    this.PS(pinds(pi)) = flag;
    this.Pinc(pinds(pi)) = true;
end
% pinds = find(this.PS);
% group = this.PS(pinds);
% layer = this.PL(pinds);
% phase = this.PP(pinds);
% pinc = this.Pinc(pinds);
% sz = [length(pinds) 5];
% varTypes = {'uint16','uint16','uint16','double','logical'};
% varNames = {'Group','Pind','Layer','Phase','Include'};
% this.phaseTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
% tab = {group,pinds,layer,phase,pinc};
% this.phaseTable = {group,pinds,layer,phase,pinc};

end



