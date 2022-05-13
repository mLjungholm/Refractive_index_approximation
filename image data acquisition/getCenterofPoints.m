% Find the mean position of a numer of points within a set span.
% pointList - cell array containing points & coordinates
% span [x1,x2] - start and end of span to be checked


% This function is unused

function centerpoint = getCenterofPoints(pointList,span)

lineNums = length(pointList);
centerpoint = zeros(lineNums,1);

for ind = 1:lineNums
    
    test_line = pointList{ind}(:,2);
    [Inds,~] = find(test_line > span(2));
    test_line(Inds) = nan;
    [Inds,~] = find(test_line < span(1));
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
end