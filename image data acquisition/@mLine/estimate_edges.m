function estimate_edges(this,options)
arguments
    this
    options {mustBeMember(options,{'left','right','all'})} = 'all'
end

gradPoints = zeros(size(this.points));
for layerInd = 1:this.imNums
    grad_temp = gradient(this.gaussPoints(:,layerInd)); % First gradient
    grad_temp = abs(gradient(grad_temp)); % Abs value of second gradient
    gradPoints(:,layerInd) = grad_temp;
end

switch options
    case 'left'
        this.plotData(interpolation = 'gauss')
        fprintf('Set start of span /n')
        [xs,~] = getXfromPlot;
        fprintf('Set end of span /n')
        [xe,~] = getXfromPlot;
        
        leftEdgeList = zeros(this.imNums,1);
        for layerInd = 1:this.imNums
            filtered_points = gradPoints(:,layerInd);
            filtered_points(this.d < xs) = 0;
            filtered_points(this.d > xe) = 0;
            [~, edgeInd] = max(filtered_points);
            leftEdgeList(layerInd) = edgeInd;
        end
        this.leftEdge = mean(leftEdgeList);
        distToEdge = (this.d-this.leftEdge);
        this.leftEdgeIndex = find(distToEdge > 0,1)-1;
        
    case 'right'
        this.plotData(interpolation = 'gauss')
        fprintf('Set start of span /n')
        [xs,~] = getXfromPlot;
        fprintf('Set end of span /n')
        [xe,~] = getXfromPlot;
        
        rightEdgeList = zeros(this.imNums,1);
        for layerInd = 1:this.imNums
            filtered_points = gradPoints(:,layerInd);
            filtered_points(this.d < xs) = 0;
            filtered_points(this.d > xe) = 0;
            [~, edgeInd] = max(filtered_points);
            rightEdgeList(layerInd) = edgeInd;
        end
        this.rightEdge = mean(rightEdgeList);
        distToEdge = (this.d-right.leftEdge);
        this.rightEdgeIndex = find(distToEdge > 0,1);
    
    case 'all'
        this.estimate_edges('left')
        this.estimate_edges('right')
end


end