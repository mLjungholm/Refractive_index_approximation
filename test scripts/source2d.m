% 2D source of rays for GRIN ray trace
classdef source2d < handle
    properties
        P = [];     % Current point vector
        V = [];     % Current direction vector
        P0 = [];    % Starting vector
        PT = {};    % Path vector
        ended = []; % Termination vector
        nRays = 0;  % Number of rays
        sWidth = 0; % Width of source
        nStart = 1; % Refractive index of start medium
        totalPath = [];
        OP = [];    % Optical path length
        phase = [];
        backTrace = [];
    end
    methods
        function this = source2d(startP,startV, nRays, sWidth,refractiveIndex)
            this.nRays = nRays;
            this.sWidth = sWidth;
            this.V = ones(nRays,2).*startV;
            this.ended = zeros(nRays,1);
            this.PT = cell(nRays,1);
            this.nStart = refractiveIndex;
            this.OP = zeros(nRays,1);
            this.phase = zeros(nRays,1);
            rg = linspace(-sWidth/2,sWidth/2,nRays);
            this.P = ones(nRays,2).*startP;
            this.P(:,2) = this.P(:,2)+rg';          
            theta = -acos(startV(1)./sqrt(startV(1)^2+startV(2)^2));
            this.backTrace = zeros(nRays,2);
            if startV(2) < 0
                theta = -theta;
            end
            if theta ~= 0
                this.P = this.P - startP;
                rot_theta = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                this.P = this.P*rot_theta;
                this.P = this.P + startP;
                this.P0 = this.P;
            end
        end
        
        function plotP(this,figureNr)
            figure(figureNr)
            hold on
            axis equal
            xlabel('X - axis')
            ylabel('Y - axis')
            plot(this.P(:,1),this.P(:,2),'.')
            quiver(this.P(:,1),this.P(:,2),this.V(:,1),this.V(:,2),'color','b')
            grid on
        end
        
        function plotTrace(this, figureNr,color)
            figure(figureNr)
            hold on; axis equal
            for rayInd = 1:this.nRays
                if ~isempty(this.PT{rayInd})
                    plot(this.PT{rayInd}(:,1),this.PT{rayInd}(:,2),color)
                end
            end
        end
        
        function totalPath = getTotalPath(this)
            totalPath = zeros(this.nRays,1);
            totalPath(:) = nan;
            for rayInd = 1:this.nRays
                if ~isempty(this.PT{rayInd})
                    tl = 0;
                    p1 = this.PT{rayInd}(1,:);
                    pNum = 2;
                    while ~isnan(this.PT{rayInd}(pNum,1))
                        p2 = this.PT{rayInd}(pNum,:);
                        p21 = p2-p1;
                        l = sqrt(p21(1)^2 + p21(2)^2);
                        tl = tl + l;
                        p1 = p2;
                        pNum = pNum + 1;
                    end
                    totalPath(rayInd) = tl;  
                end
                this.totalPath = totalPath;
            end
        end
        
        function pathHist(this,figureNr)
            figure(figureNr)
            histogram(this.totalPath,20,'binLimits',[min(this.totalPath),max(this.totalPath)])
        end
        
        function plotBacktrack(this,stopLine,figureNr)
            figure(figureNr)
            hold on, axis equal
            for rayInd = 1:this.nRays
                t = (stopLine - this.P(rayInd,2))/this.V(rayInd,2);
                px = this.P(rayInd,1) + t*this.V(rayInd,1);
                plot([this.P(rayInd,1);px],[this.P(rayInd,2),stopLine],'--b')
            end
        end
        
        function [backTrace] = getBacktrace(this,backLine)
            for rayInd = 1:this.nRays
                t = (backLine - this.P(rayInd,2))/this.V(rayInd,2);
                px = this.P(rayInd,1) + t*this.V(rayInd,1);
                this.backTrace(rayInd,:) = [px,backLine];
            end
            backTrace = this.backTrace;
        end
        
        function plotOP(this,figureNr)
            figure(figureNr)
            hold on;
            scatter(this.P0(:,1),this.phase)
        end
        
        function phaseShift = getPhaseShiftValues(this,lambda,r,xRange)
            phaseShift = this.phase.*r;
            phasePos = this.backTrace(:,1);
            phasePos = phasePos(phaseShift > 0);
            phaseShift = phaseShift(phaseShift > 0);
            phaseShift = phaseShift(phasePos > xRange(1));
            phasePos = phasePos(phasePos > xRange(1));
            phasePos = phasePos.*r;            
%             phaseShift = phaseShift(phasePos > xRange(1));
%             phasePos = phasePos(phasePos > xRange(1));
%             phaseShift = phaseShift(phasePos < xRange(2));
%             phasePos = phasePos(phasePos < xRange(2));
%             phasePos = phasePos.*r;
            phaseFunc = @(x) (cos(x*2*pi/lambda + pi) + 1)/2;
            phaseShift = [phaseFunc(phaseShift),phasePos];
            phaseShift = sortrows(phaseShift,2);
        end
    end
end



