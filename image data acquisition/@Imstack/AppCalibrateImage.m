function roiLine = AppCalibrateImage(this,imhandle)
roiLine = drawline;          % Calls draw function
tempLine = roiLine.Position; % Store the coordinates form the temporary roi object

linePixelDist = norm(tempLine(2,:)-tempLine(1,:));
this.pixelSize = (knownLineDist*mScale)/linePixelDist;
end