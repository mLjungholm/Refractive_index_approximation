%% Smooths the datapoints with a Gaussian and a Savitzky-Golay filter
function smoothData(this)
this.gaussPoints = zeros(size(this.points));
this.sgolayPoints = zeros(size(this.points));

for ind = 1:this.imNums
    l_movingmadian = smoothdata(this.points(:,ind),'movmedian',5);
    l_gaussian = smoothdata(this.points(:,ind),'gaussian',50);
    l_sgolay = smooth(l_movingmadian,40,'sgolay');
    this.gaussPoints(:,ind) = l_gaussian;
    this.sgolayPoints(:,ind) = l_sgolay;
end
end