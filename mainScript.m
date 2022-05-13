%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main script for estimating the refractive index of curved interference %
% microscopy objects                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step:1   Create and import image stack
lens1 = Imstack;
lens1.importStack;
lens1.dispRawImage(1)


%% step:2   Crop stack around lens
lens1.cropStack
lens1.dispImage(1)

%% step:3   Import Calibration image
lens1.importCalibrationImage

%% Step:4  Set variables
lambda = 550*10^-9;
n0 = 1.36;
lens1.lambda = lambda;

%% Step:5 Create any support lines if needed.
lens1.createSupportLine;

%% Step 6-? samples the image and estimates the refractive index. Repeat as many times as nessesary

%% step:6   Create sample line
lens1.drawLine

%% step:7 Set variables for the sample line (CHange id to the cvurent line)
id = lens1.lastLineId;
lens1.mLines{id}.lambda = lambda;
lens1.mLines{id}.pixelSize = lens1.pixelSize;
lens1.mLines{id}.n0 = n0;

%% step:8 Visualize sample data from the line

% lens1.mLines{id}.plotData(layerIndex = 1)
lens1.mLines{id}.plotData()
%% step:9 Determine what side to use (left or right)
% If (right) then the data will be fliped.

msg = 'Use left or right side :';
prompt = strcat(string(msg)," ",'[L/R]:');
awn = input(prompt,'s');

if isequal(awn,'r') || isequal(awn,'R')
    lens1.mLines{id}.points = flipud(lens1.mLines{id}.points);
    lens1.mLines{id}.gaussPoints = flipud(lens1.mLines{id}.gaussPoints);
    lens1.mLines{id}.sgolayPoints = flipud(lens1.mLines{id}.sgolayPoints);
    lens1.mLines{id}.d = lens1.mLines{id}.d(end) - lens1.mLines{id}.d;
    lens1.mLines{id}.d = flipud(lens1.mLines{id}.d);
    lens1.mLines{id}.dataFliped = true;
end

%% step:10 Create support lines (center and edges)

lens1.mLines{id}.normalize_points;
lens1.mLines{id}.estimate_edges('left');
close()

%% step:11 Find all peaks
lens1.mLines{id}.findPeaks

%% step 12: Estimate the center
% lens1.mLines{id}.plotData(interpolation = 'gauss')
% lens1.mLines{id}.plotSupportLines
% lens1.mLines{id}.findPeaks
lens1.mLines{id}.estimate_center
close()

%%
lens1.mLines{id}.plotPeaks(interpolation = 'gauss')
%%
lens1.mLines{id}.plotData(interpolation = 'gauss')

%%
threshold = 50;
minPeakNumber = 3;
lens1.mLines{id}.tracePeaks(threshold, minPeakNumber)
lens1.mLines{id}.centerZone = 40;
lens1.mLines{id}.find_minmax_slice
%%
lens1.mLines{id}.cycleData

%%

lens1.mLines{id}.estimate_phase_shift
%%
lens1.mLines{id}.plotPhase
