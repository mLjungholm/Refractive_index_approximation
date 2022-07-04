function clearTracedPeaks(this)
this.PD = []; % Point distance
this.PL = []; % Point Layer
this.PS = []; % Point series
this.PP = []; % Point Phase
this.PV = []; % peak value
this.L = {}; % Each cell in L(layer) = [points];
this.S = {}; % Series
this.Pinc = []; % Include point in calulations & plots
this.seriesNums = 0;
end