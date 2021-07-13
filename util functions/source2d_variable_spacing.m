% 2D source of rays for GRIN ray trace
classdef source2d_variable_spacing < handle
    properties
        P = [];     % Current point vector
        V = [];     % Current direction vector
        V0 = [];
        P0 = [];    % Starting vector
        PT = {};    % Path vector
        ended = []; % Termination vector
        nRays = 0;  % Number of rays
        nStart = 1; % Refractive index of start medium
        totalPath = [];
        OP = [];    % Optical path length
        phase = [];
        backTrace = [];
        phasePath = {};
        minR = [];
        diviation = {};
        stepError = {};
        nPath = {};
    end
    methods
        function this = source2d_variable_spacing(startP,startV,refractiveIndex)
            this.nRays = size(startP,1);
            this.P0 = startP;
            this.P = startP;
            this.V = ones(this.nRays,2).*startV;
            this.V0 = this.V;
            this.ended = zeros(this.nRays,1);
            this.PT = cell(this.nRays,1);
            this.nStart = refractiveIndex;
            this.OP = zeros(this.nRays,1);
            this.phasePath = zeros(this.nRays,1);
            this.phase = zeros(this.nRays,1);
            this.totalPath = zeros(this.nRays,1);
            this.backTrace = zeros(this.nRays,2);
            this.diviation = cell(this.nRays,1);
            this.stepError = cell(this.nRays,1);
            this.nPath = cell(this.nRays,1);
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
        
        function resetSource(this)
            this.P = this.P0;
            this.V = this.V0;
            this.ended = zeros(this.nRays,1);
            this.PT = cell(this.nRays,1);
            this.OP = zeros(this.nRays,1);
            this.phasePath = zeros(this.nRays,1);
        end
        
        function [backTrace] = getBacktrace(this,backLine)
            for rayInd = 1:this.nRays
                t = (backLine - this.P(rayInd,2))/this.V(rayInd,2);
                px = this.P(rayInd,1) + t*this.V(rayInd,1);
                this.backTrace(rayInd,:) = [px,backLine];
            end
            backTrace = this.backTrace;
        end
        function plotBacktrack(this,figureNr)
            figure(figureNr)
            hold on, axis equal
            for rayInd = 1:this.nRays
%                 t = (stopLine - this.P(rayInd,2))/this.V(rayInd,2);
%                 px = this.P(rayInd,1) + t*this.V(rayInd,1);
                plot([this.P(rayInd,1);this.backTrace(rayInd,1)],[this.P(rayInd,2);this.backTrace(rayInd,2)],'--b')
            end
        end
        function plotRay(this,rayInd,color,figureNr)
            figure(figureNr)
            plot(this.PT{rayInd}(:,1),this.PT{rayInd}(:,2),color)
        end
        function px = plotRay_back(this,rayInd,backLine,color,figureNr)
            figure(figureNr)
            t = (backLine - this.P(rayInd,2))/this.V(rayInd,2);
            px = this.P(rayInd,1) + t*this.V(rayInd,1);
            plot([this.P(rayInd,1); px], [this.P(rayInd,2); backLine],color);
        end
        function [phase,totalPath,v,p] = getEndVals_single(this,rayInd,printval,scale)
            phase = this.phase(rayInd);
            totalPath = this.totalPath(rayInd);
            v = this.V(rayInd,:);
            p = this.P(rayInd,:);
            if isequal(printval,'print')
                fprintf('V = [%.6f, %.6f], P = [%.6f, %.6f], Total phase = %.4f, Total path = %.4f \n',v(1),v(2),p(1),p(2),phase/scale,totalPath/scale)
            end
        end
    end
end




