function calibrateImage(this)
if isempty(this.calibIm)
    fprintf('Error: there is no calibration image. Use Imstack.importCalibrationImage() \n')
    return
end

fprintf('Draw a line in the image for a known distance. \n Press any key when done \n')
figure('Name','Draw distance','NumberTitle','off')
imshow(this.calibIm)
set(gcf,'position',this.figSize)
axis equal

roiLine = drawline;          % Calls draw function
pause                        % Waits for key press confirmation
tempLine = roiLine.Position; % Store the coordinates form the temporary roi object
close();

flag = 1;
while flag
msgIn = input('What is the scale used (nm, um, mm...) : ','s');
switch msgIn
    case 'nm'
        mScale = 10^-9;
        flag = 0;
    case 'um'
        mScale = 10^-6;
        flag = 0;
    case 'mm'
        mScale = 10^-4;
        flag = 0;
    case 'cm'
        mScale = 10^-2;
        flag = 0;
    otherwise
        mScale = nan;
        fprintf('Error: %s is not a valid scale /n',msgIn)
end
knownLineDist = input('Input the distance for the calibration line: ');
linePixelDist = norm(tempLine(2,:)-tempLine(1,:));
this.pixelSize = (knownLineDist*mScale)/linePixelDist;
end