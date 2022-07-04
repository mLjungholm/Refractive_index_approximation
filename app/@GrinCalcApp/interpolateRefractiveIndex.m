function interpolateRefractiveIndex(app)
simLines = cell2mat(app.SamplingLineTable.Data(:,3)) & cell2mat(app.SamplingLineTable.Data(:,4));
simLinesN = nnz(simLines);
if simLinesN < 2
    return
end
% inices of the sampling lines
lineInds = find(smiLines);

% Number of pints in the lines
% lineLength = length(app.lens.mLines{lineInds(1)}.gradientD);

% Create an array with the x,y coordinates and n-values for each line
pointList = cell(simLinesN,1);
for lind = 1:simLinesN
    % Generate x,y coordinates
    xy = app.lens.mLines{lineInds(lind)}.npoints2imcoords(); % [x,y} for the sampling points
    pointList{lind} = [xy, app.lens.mLines{lineInds(lind)}.gradientD];
end
pointList = cell2mat(pointList);
% Find the convex hull of the lines
F = scatteredInterpolant(pointList(:,1),pointList(:,2),pointList(:,3));
k = convhull(pointList(:,1),pointList(:,2));


xgrid = min(pointList(:,1)):1:max(pointList(:,1));
ygrid = min(pointList(:,2)):1:max(pointList(:,2));
overlay_coords = cell(length(xgrid),length(ygrid));

for xind = 1:length(xgrid)
    for yind = 1:length(ygrid)
        overlay_coords{xind,yind} = [xgrid(xind),ygrid(yind)];
    end
end
overlay_coords = reshape(overlay_coords,[],1);
overlay_coords = cell2mat(overlay_coords);

[in,on] = inpolygon(overlay_coords(:,1),overlay_coords(:,2),k(:,1),k(:,2));
overlay_n = F(overlay_coords(:,1),overlay_coords(:,2));
overlay_n(~(in | on)) = nan;

end