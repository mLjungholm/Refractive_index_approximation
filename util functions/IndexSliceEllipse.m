% In use
classdef IndexSliceEllipse < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id = 'Undefined'
        gradient
        n
        nh
        r
        rMinor
        n0
        nShells
        mPhaseRadius
        mPhaseShift
        analyticals
        analyticsNames = ["nApprox" "nHigh" "nLow" "totalPhase"...
           "backPhase" "innerPhase" "outerPhase" "accuarcy" "itterations"];
    end
    
    methods
        function this = IndexSliceEllipse(mPhaseRadius,mPhaseShift,n0,rMinor)
            this.mPhaseRadius = mPhaseRadius;
            this.mPhaseShift = mPhaseShift;
            this.n0 = n0;
            this.nShells = length(mPhaseRadius)-1;
            this.r = mPhaseRadius;
            this.rMinor = rMinor;
            this.gradient = zeros(this.nShells,1).*nan;
            this.n = zeros(this.nShells + 1,1).*nan;
            this.n(1) = n0;
            this.nh = zeros(this.nShells + 1,1).*nan;
            this.nh(1) = n0;
            this.analyticals = cell(this.nShells,length(this.analyticsNames));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                      Utility functions                          %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %-----------------------------------------------------------------%
        % Get the refractive index for a point with the gradient in shell i
        function n = getN(this, shell_index, p)
            if isnan(this.n(shell_index))
                fprintf('Error in IndexSlice "%s" \n: n(%u) is not defiened\n',this.id, shell_index) 
                return
            elseif isnan(this.gradient(shell_index))
                fprintf('Error in IndexSlice "%s" \n: gradient(%u) is not defiened\n',this.id, shell_index)
                return
            end
            theta = acos(p(1)/norm(p));
            r0 = this.r(1)^2*this.rMinor/sqrt((this.r(1)*sin(theta))^2 + (this.r(1)*this.rMinor*cos(theta))^2);
            rp = norm(p)*this.r(1)/r0;
            n = this.gradient(shell_index) * (rp - this.r(shell_index)) + this.n(shell_index);            
        end
        
        
        %-----------------------------------------------------------------%                
        % Calculate the gradient in a shell or all shells
        function calculateGradient(this,shell)
            if isequal(shell,'all') 
                for i = 1:this.nShells
                    k = (this.n(i)-this.n(i+1))/(this.r(i)-this.r(i+1));
                    this.gradient(i) = k;
                end
            else
                if isnan(this.n(shell)) || isnan(this.n(shell+1))
                    fprintf('Error in IndexSlice "%s" \n: n(%u) is not defiened\n',this.id, shell)
                    return
                end
                k = (this.n(shell)-this.n(shell+1))/(this.r(shell)-this.r(shell+1));
                this.gradient(shell) = k;
            end
        end
        
        %-----------------------------------------------------------------%
        % Calculate the phase-shift at a given distance r
        function phaseShift = getPhaseShift(this,p)
            % Returns the index of mesured phase positions with the smalest
            % radius larger than "r"
            rp = norm(p);
            theta = acos(p(1)/rp);
            r0 = this.r(1)^2*this.rMinor/sqrt((this.r(1)*sin(theta))^2 + (this.r(1)*this.rMinor*cos(theta))^2);
            rp = rp*this.r(1)/r0;
            pind = find((this.mPhaseRadius-rp) > 0, 1,'last'); 
            
            if isempty(pind) % No radius larger than "r"
                pind = 1;
                
            elseif pind == length(this.mPhaseRadius) % All radius larger than "r"
                pind = length(this.mPhaseRadius) - 1;
            end
            
            p0 = this.mPhaseRadius(pind);
            p1 = this.mPhaseRadius(pind+1);
            
            shift0 = this.mPhaseShift(pind);
            shift1 = this.mPhaseShift(pind+1);
            
            phaseShift = (shift1-shift0)/(p1-p0)*(rp-p0) + shift0;            
            if phaseShift < 0
                phaseShift = 0;
            end
        end
        
        %-----------------------------------------------------------------%
        % Calculate the current radius fraction of R-major axis
        function rr = getRelativeRadius(this,p)
            rp = norm(p);
            theta = acos(p(1)/rp);
            r0 = this.r(1)^2*this.rMinor/sqrt((this.r(1)*sin(theta))^2 + (this.r(1)*this.rMinor*cos(theta))^2);
            rr = rp*this.r(1)/r0;
        end       
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                   Visualization functions                       %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %-----------------------------------------------------------------%
        % Plot the refractive index for a shell or all
        function plotN(this,varargin)
            if isempty(varargin)                
                scatter(this.r(2:end),this.n(2:end),'b*')
            else
                scatter(this.r(varargin{1}+1),this.n(varargin{1}+1),'b*')
            end
        end
        
        
        %-----------------------------------------------------------------%
        % Plot the refractive index of all shells as line segments
        function plotNCurve(this)
            for i = 1:this.nShells
                x = [this.r(i); this.r(i+1)];
                y = [this.n(i); this.n(i+1)];
                plot(x,y,'r')
            end
        end
        
        %-----------------------------------------------------------------%
        % Plot the refractive index of all shells as line segments
        function plotGradient(this,varargin)
            if isempty(varargin)
                for i = 1:this.nShells
                    x = [this.r(i); this.r(i+1)];
                    y = [this.n(i); this.getN(i,this.r(i+1))];
                    plot(x,y,'r')
                end
            else
                for i = 1:length(varargin{1})
                    ind = varargin{1}(i);
                    x = [this.r(ind); this.r(ind+1)];
                    y = [this.n(ind); this.getN(ind,this.r(ind+1))];
                    plot(x,y,'r')
                end
            end
        end
        
        %-----------------------------------------------------------------%
        % Plot the measured phase profile
        function plotPhase(this)
            plot(this.mPhaseRadius, this.mPhaseShift,'b')
        end
        
        %-----------------------------------------------------------------%
        % Plot refractive index analytics
        function plotAnalyticsN(this,rayInd)
            x = 1:length(this.analyticals{rayInd,1});
            plot(x,this.analyticals{rayInd,1},x,this.analyticals{rayInd,2},x,this.analyticals{rayInd,3})
        end
    end
end
