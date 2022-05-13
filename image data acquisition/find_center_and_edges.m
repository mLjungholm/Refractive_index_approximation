% Script for finding the center of the double sided curve and the edges of
% the lens.

% The function uses the lineVals and LinePeaks arrays as input.


%% Find center
% min, max for distance coordinate
center_intervall = [240 300];

centerpoint = zeros(lineNums,1);
for ind = 1:lineNums
    
    test_line = linePeaks{ind}(:,2);
    [Inds,~] = find(test_line > center_intervall(2));
    test_line(Inds) = nan;
    [Inds,~] = find(test_line < center_intervall(1));
    test_line(Inds) = nan;
    test_coord = test_line(~isnan(test_line));
    
    if isempty(test_coord)
        test_coord = nan;
    elseif length(test_coord) > 1
        test_coord = sum(test_coord)/length(test_coord);
    end
    centerpoint(ind) = test_coord;
end

centerpoint = centerpoint(~isnan(centerpoint));
centerpoint = sum(centerpoint)/length(centerpoint);


%% Find left edge

