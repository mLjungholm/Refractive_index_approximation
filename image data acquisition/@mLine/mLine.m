%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         Sampling line class                             %
%                                                                         %
% Class that holds the information of a sampling line. It contains all    %
% sampling data for the line, functionf for modifying the data aswell as  %
% information on peak locations and phase shift data.                     %
%                                                                         %
% Version 1.1                                                             %
% By: Mikael Ljungholm, 2020/03                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef mLine < handle
    
    properties
        lineInd = nan;      % Line id
        lineCoord = [];     % [x1,y1; x2, y2] for the line
        points = [];        % sampled data points [point_index, layer_index]
        normalized_points = 0; % Flag for identifying if the points have been normalized
        coords = [];        % [x,y] for the points
        d = [];             % Distance from start for the points
        pixelSize = nan;    % True size of pixels
        lambda = nan;       % Wavelength
        n0 = nan;           % Initial refractive index.
        refractiveGradient = []; % Calculated refractive index gradient.
        gradientD = [];
        pointNums = nan;    % Number of sampling points
        imNums = nan;       % Number of images in the sampling stack
        testSlice = nan;
        
        gaussPks = [];          % peak coordinate for the gaussian
        gaussPoints = [];       % Gaussian point coordinates
        sgolayPks = [];         % peak coordinates for the Savitzky-Golay {pointInd,layer}
        sgolayPoints = [];      % Sgolay point coordinates
        pks = [];               % summarized peaks
        pksNums = nan;          % Number of summarized peaks
        phasePoints = [];
        phase_func = nan;
        fittype = 'none';
        phaseFitChangeFlag = false;
        
        centerLine = nan;       % Coordinate for the center of the lens
        centerLineIndex = nan;
        leftEdge = nan;         % Edge of lens in relation to plot
        leftEdgeIndex = nan;
        rightEdge = nan;        % -||-
        rightEdgeIndex = nan;
        leftPhaseMax = nan;    % Index for the layer used as reference point for the phase shift
        leftPhaseMin = nan;
        
%         sgolayZone = 0;
%         centerVal = [];
        
        centerZone = nan;
        edgeZone = nan;
        sgolayEdge = nan;
%         sgolayCenter = nan;
        
        %         Variables from the trace peaks function
        PD = [] % Point distance
        PL = [] % Point Layer
        PS = [] % Point series
        PP = [] % Point Phase
        PV = [] % peak value
        L = {} % Each cell in L(layer) = [points];
        S = {} % Series
        Pinc = [] % Include point in calulations & plots
        seriesNums = 0;
        layerPhaseStep = nan;
        edgeFit = nan;
%         phaseTable = table()
        
        validS = [];
        validP = [];
        validPS = [];
        
        dataFliped = false;
    end
    
    methods
        % Instantiation function
        function this = mLine(coords, lineInd)
            this.lineCoord = coords;
            this.lineInd = lineInd;
%             sz = [1 5];
%             varTypes = {'uint8','uint8','uint8','double','logical'};
%             this.phaseTable = table('Size',sz,'VariableTypes',varTypes);
        end
        
        % Function block
        plotData(this,options) %ok%    % Plot sampled data
        plotPeaks(this, options)    % Plot the peaks for one layer
        plotPhase(this)
        draw_line_on_image(this) % Shift this to the imstack class
        calculate_point_distance(this)
        estimate_center(this,type)
        
        smoothData(this) % ok % Smooth sampled data (savitzky-Golay & Gaussian)
        setEdgeManualy(this, side) % ok %
        
        findPeaks(this) %ok%    % Try to find the peaks from the data
        plotSupportLines(this) % Plot center and edge lines if any        
        
        cycleData(this)   %ok%
        tracePeaks(this,threshold, minPeakNumber) % in development
        estimate_phase_shift(this)
        normalize_points(this) % Normalizes the values of each distance point for the total interwall for that point.
        reset_line(this,options) % options = {'all','peaks','edges','phase'} default = 'all'
        estimate_edges(this,options) % options = {'left','right','all'} default = 'all'
        removeOutlier(this)
        fit_phase_shift_profile(this)
        validation_trace(this)
        calculate_refractive_index(this)
        grin = create_grin(this, stepNums)
        find_minmax_slice(this)
        clearTracedPeaks(this)
        xy = npoints2imcoords(this) % Converts the refractive vakues to the x,y coordinates from the sampling image

        AppPlotData(this,imhandle,interpolation,layerIndex)
        AppPlotPeaks(this,imhandle,interpolation,layerIndex)
        AppFlipData(this)
        AppPlotPhaseGradient(this, i)
        AppEstimatePhaseShift(this)
        flag = AppTracePeaks(this,layer,span,threshold, min_nums_in_serie)
        AppCalculateRefractiveIndex(this)
        AppFitPhaseShiftProfile(this,fittype,lowRes)
        AppPlotRefractiveGradient(imhandle)
        % Test functions
        function plotTracedPeaks(this)
            figure(1)
            grid on
            hold on
            for ind = 1:length(this.pks)
                plot(this.pks{ind},ones(length(this.pks{ind})).*ind,'o')
            end
        end
    end
end









