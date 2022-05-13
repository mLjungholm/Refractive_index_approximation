% Function that tries to trace all the peaks thorugh all layers.
% 
% ! In development. Right now it only works for the left side and gaussian
% smoothing
function tracePeaks(this, threshold, minPeakNumber)

% threshold = 10; % max distance between two points
min_nums_in_serie = minPeakNumber; % minimum number of consecutive points in a series

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
numsS = 0; % number of series

% P = cell(1,numsP); % Each cell in P(points) = [dist, layer, series, phase];
% P = cellfun(@(x) zeros(4,1),P,'UniformOutput',false);
this.PD = zeros(numsP,1);
this.PL = zeros(numsP,1);
this.PS = zeros(numsP,1);
this.PP = zeros(numsP,1);
this.PV = zeros(numsP,1); % peak value
this.L = cell(1,numsL); % Each cell in this.L(layer) = [points]; 
this.S = {}; % Repository for all peak series this.S(series) = [points];
checkedPeaks = zeros(numsP,1);


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


for startP = 1:numsP
    if checkedPeaks(startP)
        continue
    end
    traceList = zeros(numsL,1); % Create empty list (one possible slot for each layer)
    traceList(this.PL(startP)) = startP; % Add the point index to the corresponding layer
    nums_in_serie = 1;  % Keep track of the numer of points in the series
    checkedPeaks(startP) = 1;   % Log the current point as checked
    start_layer = this.PL(startP);   % What layer does the trace start at
%     current_layer = start_layer; % What is the current layer to be tested.
    
    % If the image sequence aw perfect then you should be able to loop
    % around the last image and contiue on the first. In this case it is
    % not so. So the code will end att the last or first layer.
    d0 = this.PD(startP); % Set the first point distance
    for layer_ind = start_layer:numsL-1
        test_layer = layer_ind + 1; % Test layer
        points_in_layer = this.L{test_layer};    % all point indices in that layer
        dist_in_layer = this.PD(points_in_layer);    % Corresponding distance for each point
        [deltaP, closest_ind] = min(abs(d0 - dist_in_layer));   % Find the point that is closest
%         fprintf('Delta d = %3.2f \n',deltaP)
        closest_point = points_in_layer(closest_ind); % Get the global index for that point  
        if deltaP < threshold && ~checkedPeaks(closest_point) && this.PD(closest_point) >= d0  % Check if the distance is below the threshold and if the point has been checked.
            traceList(test_layer) = closest_point; % Add point to list
            nums_in_serie = nums_in_serie + 1; % Increase count
            checkedPeaks(closest_point) = true; % The point has now been checked.
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
        if deltaP < threshold && ~checkedPeaks(closest_point) && this.PD(closest_point) <= d0  % Check if the distance is below the threshold and if the point has been checked.
            traceList(test_layer) = closest_point; % Add point to list
            nums_in_serie = nums_in_serie + 1; % Increase count
            checkedPeaks(closest_point) = true; % The point has now been checked.
            d0 = this.PD(closest_point);
        else
            break
        end
    end
    if nums_in_serie >= min_nums_in_serie
        numsS = numsS + 1;
        this.S{numsS} = traceList(traceList ~= 0);
    end
end

% Add what series each point belongs to
for s_ind = 1:numsS
    p_inds = this.S{s_ind};
    this.PS(p_inds) = s_ind;
end
% Correct the distances to the left edge
this.PD = this.PD - this.leftEdge;


% This code is for testploting a single series of peaks

% test_ind = 1;
% close all
% 
% figure(1)
% hold on
% grid minor
% plot(this.PD(this.S{test_ind}),ones(length(this.S{test_ind}),1).*test_ind,'.')
% 
% figure(2)
% hold on
% grid minor
% for i = 1:length(this.S{test_ind})
%     point_ind = this.S{test_ind}(i);
%     lineInd = this.PL(point_ind);
%     %     plot_pks(this,line_ind,point_ind,this.PD,this.PV)
%     plot(this.d,this.gaussPoints(:,lineInd),'color',[0 0.4470 0.7410].*0.8)
%     if point_ind ~= 0
%         plot(this.PD(point_ind),this.PV(point_ind),'o','color',[0 0.4470 0.7410].*0.4)
%     else
%         plot(this.pks{lineInd}(:,2),this.pks{lineInd}(:,1),'o','color',[0 0.4470 0.7410].*0.4)
%     end
%     pause
% end
end
