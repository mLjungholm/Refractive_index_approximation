function importCalibrationImage(this)
if ~isempty(this.calibIm)
    fprintf('A calibration image already exsists \n')
    flag = yesno('Do you want to overwrite the calibration image?');
    if ~flag
        return
    end
end

fprintf('Select a calibration image \n')
[filename, path] = uigetfile('*');
filePath = strcat(path,filename);
I = imread(filePath);
this.calibIm = I;


flag = yesno('Do you want to calibrate the image?');
if ~flag
    return
end
this.calibrateImage
end