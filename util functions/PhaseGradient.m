classdef PhaseGradient < handle
    properties
        measuredPhase = [];
        phasePosition = [];
        measuredPoints = 0;
        phaseGradient = [];
        nShells = 0;
        n = [];
        shellR = [];
        r = 0;
        n0 = 0;
    end
    
    methods
        function this = PhaseGradient(phasePosition,measuredPhase,n0,r)
            this.measuredPhase = measuredPhase;
            this.phasePosition = phasePosition;
            this.measuredPoints = length(measuredPhase);
            this.shellR = this.phasePosition(this.phasePosition < r);
            this.n  = zeros(length(this.shellR),1);
            this.n0 = n0;
            this.r = r;
        end
        
        function createGradient(this)
            if ~nnz(this.n)
                disp('Error: no registered refracitve index')
                return
            end
            dr = this.shellR(1) - this.r;
            dn = this.n(1) - this.n0;
            this.phaseGradient(1) = dn/dr;
            for shellInd = 1:this.nShells-1
                dr = this.shellR(shellInd + 1) - this.shellR(shellInd);
                dn = this.n(shellInd + 1) - this.n(shellInd);
                this.phaseGradient(shellInd+1) = dn/dr;
            end
        end
        
        function plotPhase(this,figureNr)
            figure(figureNr)
            hold on; axis equal
            plot(this.measuredPhase,this.phasePosition)
        end
        
        function plotN(this,figureNr)
            if ~nnz(this.n)
                disp('Error: no registered refracitve index')
                return
            end
            figure(figureNr)
            hold on; axis equal
            plot(this.n, this.shellR)
        end
        
    end
end