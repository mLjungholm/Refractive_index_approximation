%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         Image analysis class                            %
%                                                                         %
% Main analysis class. THe Imstack object holds information on the        %
% imported images and a container for all defined sampling lines as well  %
% as visualization functions and analysis functions.                      %  
%                                                                         %
% Version 1.0                                                             %
% By: Mikael Ljungholm, 2020/02                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef  Imstack < handle
    
    properties
        imStackRaw = [];    % Stack of imported images (unmodifyed)
        sizeRaw = [];       % Original size of imported images
        imStack = [];       % Stack of images after croping
        imNums = nan;       % Number of images in stack
        imSize = [];        % Size of images after croping
        cropRect = [];      % Croping rectangle
        supLines = [];      % List of lines used for vizualisation support
        mLines = [];       % List of sampling lines
        lineNums = nan;     % Number of sampling lines
        viewSize = [300,100,1000,800];  % Desired size of ploted images
        figSizeRaw = nan;   % ? is this not the same as viewSize? Have I duplicated a variable?
        figSize = nan;      % Figure size after modifications to fit images.
        calibIm = [];       % Size calibration image
        calibLine = [];
        calibKnownDistance = nan;
        pixelSize = nan;    % Size of each pixel
        lambda = nan;       % wavelength of test source
        imPhaseShift = nan; % Phase shift between each image (fraction of 2pi phase shift)
%         phaseShiftMap = []; % Image containing all phase shift points and values
%         rIndexMap = [];     % Calculated refractive index map
        lastLineId = nan;
        n0 = nan;
        n_overlay = [];     % Imstack overlay with interpolated n-values
%         overlayPos = [];    % image position of the n-overlay
        includeLine = [];   % boolean if the line should be included in calculations
        simulatedLine = []; % boolean if the line has been simulated
    end
    
    methods
        testFun(this)
        drawLine(this)
        plotLines(this)
        cropStack(this,imhandle)
        createSupportLine(this)        
        importCalibrationImage(this)
        calibrateImage(this)
        calculate_refractive_index(this,lineIndex)
        
        AppCropStack(this,roi)
        AppPlotLines(this,imhandle,linetype,lineid)
        AppDeleteLine(this,linetype,id)        
        ind = AppCreateSamplingLine(this,roi)
        ind = AppCreateSupportLine(this,imhandle)
        roiLine = AppCalibrateImage(this,imhandle)
        AppInterpolateRefractiveIndex(this)
        AppUppdateLineData(this,lineData)

        %% Import a stack of images
        function importStack(this)
            [filename, path] = uigetfile('*','Multiselect','on');
            this.imNums = length(filename);
            filePath = strcat(path,filename{1});
            I = imread(filePath);
            this.sizeRaw = size(I);
            this.imStackRaw = uint8(zeros(this.sizeRaw(1),this.sizeRaw(2),this.imNums));
            for fileId = 1:this.imNums
                filePath = strcat(path,filename{fileId});
                I = imread(filePath);
                this.imStackRaw(:,:,fileId) = I;
            end
            ds = this.viewSize(3)/this.sizeRaw(2);
            if this.sizeRaw(2)*ds > 800
                ds = this.viewSize(4)/this.sizeRaw(1);
            end
            this.figSizeRaw = [this.viewSize(1) this.viewSize(2)...
                ds*this.sizeRaw(2) ds*this.sizeRaw(1)];
        end
        
        function dispImage(this, imageIndex)
            figurename = strcat('cropped image Nr: ',num2str(imageIndex));
            figure('Name',figurename,'NumberTitle','off')
            imshow(this.imStack(:,:,imageIndex))
            set(gcf,'position',this.figSize)
            axis equal
        end
        function dispRawImage(this, imageIndex,imhandle)
%             figurename = strcat('Raw image Nr: ',num2str(imageIndex));
%             figure('Name',figurename,'NumberTitle','off')
            imshow(this.imStackRaw(:,:,imageIndex),'Parent',imhandle)
%             set(gcf,'position',this.figSizeRaw)
%             axis equal
        end
        
        
        %% Delete unwanted sampling line %%
        function deleteLine(this, lineInd)
            if isnan(this.lineNums)
                fprintf('ERROR: There are no data sampling lines \n')
                return
            elseif this.lineNums > 1
                this.mLines{lineInd} = [];
                this.lineNums = this.lineNums - 1;
            else
                this.lineNums = nan;
                this.mLines = [];
            end
        end
        
        %% Delete support line
        function deleteSupportLine(this, lineInd)
            if isempty(this.supLines)
                fprintf('ERROR: There are no support sampling lines \n')
                return
            elseif isempty(this.supLines{lineInd})
                fprintf('ERROR: There is no support line with id %f2 \n',lineInd)
            else
                this.supLines{lineInd} = [];
                emptyCell = cellfun(@isempty,this.supLines);
                if ~any(emptyCell)
                    this.supLines = [];
                end
            end
        end
        
        %% Sample data function %%
        function sampleData(this,lineId)
%             this.mLines{lineId}.croped = false;
            x = this.mLines{lineId}.lineCoord(:,1);
            y = this.mLines{lineId}.lineCoord(:,2);
            [cx,cy,tempP] = improfile(this.imStack(:,:,1),x,y); % Fist line profile to get line size
            this.mLines{lineId}.coords = [cx,cy];
%             this.mLines{lineId}.coordsRaw = [cx,cy];
            dataPoints = zeros(length(tempP),this.imNums);
            dataPoints(:,1) = tempP;
            for imInd = 2:this.imNums
                [~,~,cv] = improfile(this.imStack(:,:,imInd),x,y); % Sample line for each image in stack
                dataPoints(:,imInd) = cv;
            end
%             this.mLines{lineId}.pointsRaw = dataPoints;
            this.mLines{lineId}.points = dataPoints; % Store data pointsf
%             this.mLines{lineId}.pointNumsRaw = length(tempP);
            this.mLines{lineId}.pointNums = length(tempP);
            this.mLines{lineId}.imNums = this.imNums;
%             this.mLines{lineId}.getDistance; % Call distance function

            this.mLines{lineId}.calculate_point_distance; % Call distance function
            this.mLines{lineId}.smoothData; % Call smoothing function
        end
        
        %% Extract peaks form lines
        function getPeaks(this, lineInd)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            end
            if isequal(lineInd, 'all')
                for lineInd = 1:length(this.mLines)
                    if ~isempty(this.mLines{lineInd})
                        this.mLines{lineInd}.findPeaks;
%                         for layerInd = 1:this.imNums
%                             this.mLines{lineInd}.findPeaks(layerInd)
%                         end
                    end
                end
            else
                if isempty(this.mLines{lineInd})
                    fprintf('Error: Line with index %f2 does not exist \n',lineInd)
                    return
                end
                this.mLines{lineInd}.findPeaks;
%                 for layerInd = 1:this.imNums
%                     this.mLines{lineInd}.findPeaks(layerInd)
%                 end
            end
        end
        
        %% Plot the sampled data from line
        function plotSampledData(this,id,varargin)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif isempty(this.mLines{id}) 
                fprintf('Error: Line with index %f2 does not exist \n',id)
                return
            end
            if isempty(varargin)
                this.mLines{id}.plotData;
            elseif isequal(varargin{1},'croped')
                this.mLines{id}.plotData('croped');
            elseif isequal(varargin{1},'raw')
                this.mLines{id}.plotData('raw');
            end
            
        end
        
        
        %% Plot peaks form the smoothed data
        function plotPeaks(this, lineInd, layerInd)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif isempty(this.mLines{lineInd}) 
                fprintf('Error: Line with index %f2 does not exist \n',lineInd)
                return
            end
            this.mLines{lineInd}.plotPeaks(layerInd)
        end
        
        %% Plot the peaks from a single smoothing set
        function plotDataLine_singleSet(this,lineInd,type)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif isempty(this.mLines{lineInd}) 
                fprintf('Error: Line with index %f2 does not exist \n',lineInd)
                return
            end
            this.mLines{lineInd}.plotSingleDataSet(type);
        end
                %% Plot the peaks from a single smoothing set
        function plotDataLine_singleSetRaw(this,lineInd,type)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif isempty(this.mLines{lineInd}) 
                fprintf('Error: Line with index %f2 does not exist \n',lineInd)
                return
            end
            this.mLines{lineInd}.plotSingleDataSetRaw(type);
        end
        
        %% Function for estimation of the center
        function getLineCenter(this, lineInd, type)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif isempty(this.mLines{lineInd})
                fprintf('Error: Line with index %f2 does not exist \n',lineInd)
                return
            elseif ~isequal(type,'gauss')
                if ~isequal(type,'sgolay')     
                    strType = string(type);
                    fprintf('Error: %s is not a valid type \n',strType)
                    return
                end
            end
            this.mLines{lineInd}.plotSingleDataSetRaw(type)
            this.mLines{lineInd}.estimateCenter(type);
            close
        end
        %% Setting the center line manualy
        function setCenterManualy(this, lineInd)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif length(this.mLines) < lineInd || isempty(this.mLines{lineInd})
                fprintf('Error: Line with id %f2 does not exist \n',lineInd)
                return
            end
            this.plotDataLine_singleSet(lineInd,'raw')
            [xs,~] = getXfromPlot;
            this.mLines{lineInd}.centerLine = xs;
            close
        end
        
        %% Setting the edge line manualy
        function setEdgeManualy(this, lineInd, side)
            if isempty(this.mLines)
                fprintf('Error: There are no sampling lines \n')
                return
            elseif length(this.mLines) < lineInd || isempty(this.mLines{lineInd})
                fprintf('Error: Line with id %f2 does not exist \n',lineInd)
                return
            end
            if isequal(side,'left') || isequal(side,'l')
                this.plotDataLine_singleSetRaw(lineInd,'raw')
                [xs,~] = getXfromPlot;
                this.mLines{lineInd}.leftEdge = xs;
                close
            elseif isequal(side,'right') || isequal(side,'r')
                this.plotDataLine_singleSetRaw(lineInd,'raw')
                [xs,~] = getXfromPlot;
                this.mLines{lineInd}.rightEdge = xs;
                close
            else
                sidestr = string(side);
                fprintf('Error: %s is not a valid argument \n',sidestr)
                return
            end
        end
        
        %% Set wavelength
        function setWavelength(this)
            this.lambda = input('Input wavelength used [nm] : ');
            this.lambda = this.lambda * 10^-9;
        end
        %% Set phaseshift between each image
        function setImagePhaseShift(this)
            msgIn = input('Input phaseshift between each image (in fraction of a full phase shift) \n Enter 0 if the the image stack is a full phase shift: ');
            if msgIn == 0                
                this.imPhaseShift = 2*pi/(this.imNums - 1);
            else
                this.imPhaseShift = msgIn;
            end
        end
        %% Cropp data set for a specific line
        % This requiers edge coordinates
        function cropSampledData(this, lineInd)
            if ~this.lineExists(lineInd)
                return
            end
            this.mLines{lineInd}.cropData;
        end
        
        %% Select what peaks to use
        function selectPeaks(this,lineInd, selectionMethod)
            arguments
                this
                lineInd double
                selectionMethod {mustBeMember(selectionMethod,['span','manual','gauss','sgolay'])}
            end
            if ~this.lineExists(lineInd)
                return
            end
            this.mLines{lineInd}.selectPeaks(selectionMethod)
            
        end
        
        %% Check if line exists
        function flag = lineExists(this, lineInd)
            if isempty(this.mLines{lineInd}) || lineInd > this.lineNums
                flag = false;
                fprintf('Error: The line with id %s does not exist \n',lineInd)
            else
                flag = true;
            end
        end
    end
end














