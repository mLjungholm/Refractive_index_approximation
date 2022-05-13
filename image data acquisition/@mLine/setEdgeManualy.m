
function setEdgeManualy(this, side)

if isequal(side,'left') || isequal(side,'l')
    this.plotData(interpolation = 'raw')
    [xs,~] = getXfromPlot;
    this.leftEdge = xs;
    close
elseif isequal(side,'right') || isequal(side,'r')
    this.plotData(interpolation = 'raw')
    [xs,~] = getXfromPlot;
    this.rightEdge = xs;
    close
else
    sidestr = string(side);
    fprintf('Error: %s is not a valid argument \n',sidestr)
    return
end
end