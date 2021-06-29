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
        phasePath = {};
        minR = [];
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
        function plotRay(this,rayNr,scaleRay,figureNr)
            figure(figureNr)
            hold on
            axis equal
            xlabel('X - axis')
            ylabel('Y - axis')
            plot(this.P(rayNr,1),this.P(rayNr,2),'.')
            quiver(this.P(rayNr,1),this.P(rayNr,2),this.V(rayNr,1),this.V(rayNr,2),scaleRay,'color','b')
            grid on
        end
        
        function plotTrace(this, figureNr)
            figure(figureNr)
            for rayInd = 1:this.nRays
                if ~isempty(this.PT{rayInd})
                    plot(this.PT{rayInd}(:,1),this.PT{rayInd}(:,2),'r')
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
    end
end




