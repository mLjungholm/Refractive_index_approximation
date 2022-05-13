% circular 2D GRIN surface for testing GRIN ray trace
% The surface has a radius of 1 and is centered at [0,0]
classdef GRIN2d_rotSym < handle
    
    properties
        %         grinDim = nan; % Axis dimension of GRIN
        inputGradient = nan;
        inputDist = nan;
        X = []; % x coordinates
        Y = []; % y coordinates
        %         Z = []; % z coordinates
        P = []; % index of refraction matrix
        %         gradient_type = ''; % Type of gradient
        %         data_format = '';   % data format, vector or matrix
        stepsize = [];      % grid spacing
        gridNums = nan;
        dX = [];            % nDx matrix. Dx = gradient of P(x)
        dY = [];            % nDy matrix. Dy = gradient of P(y)
        %         dZ = [];            % nDy matrix. Dz = gradient of P(z)
        nEdge = 1;          % Refractive index att edge
        %         symetryAxis = nan;  % Rotational symetry around axis (nan = no symetry)
        %         gradientType = nan; % Input data in the form of vector or matrix
    end
    
    methods
        function this = GRIN2d_rotSym(refractiveIndex, distances, gridNums)
            arguments
                refractiveIndex double
                distances double
                gridNums double = 10^3;
            end
            
            this.gridNums = gridNums;
            if isnan(refractiveIndex(end))
                dn = (refractiveIndex(end-1)-refractiveIndex(end-2))/(distances(end-1)-distances(end-2));                
                refractiveIndex(end) = refractiveIndex(end-1) + dn/2*(distances(end)-distances(end-1));
            end
            this.inputGradient = refractiveIndex;
            this.inputDist = distances;
            r0 = max(distances);
            
            %             this.stepsize = 1/gridNums;
            %             xp = (-1:this.stepsize:1);
            %             yp = (-1:this.stepsize:1);
            this.stepsize = r0/gridNums;
            xp = (-r0:this.stepsize:r0);
            yp = (-r0:this.stepsize:r0);
            [this.X,this.Y] = meshgrid(xp,yp);
            mask = sqrt(this.X.^2 + this.Y.^2);
            mask(mask>(r0+this.stepsize/2)) = nan;
            mask(~isnan(mask)) = 1;
            
            r = sqrt(this.X.^2 + this.Y.^2).*mask;
            r = reshape(r,[],1);
            
            p = interp1(distances, refractiveIndex, r, 'linear');
            p(isnan(p)) = refractiveIndex(end);
            
%             rU = unique(r);
%             rU = rU(~isnan(rU));
%             p_map = interp1(distances, refractiveIndex, rU, 'linear');
%             p = reshape(mask,[],1);
%             for i = 1:length(r)
%                 if isnan(p(i))
%                     continue
%                 else
%                     p(i) = p_map(rU == r(i));
%                 end
%             end
            
%             this.P = reshape(p,size(mask,1),[]);
            this.P = reshape(p,size(mask,1),[]);
            [px,py] = gradient(this.P,this.stepsize);
            this.dX = this.P.*px;
            this.dY = this.P.*py;
        end
        
        function plot_nIndex(this)
            figure('Name','Refractive index plot','NumberTitle','off')
            surf(this.X,this.Y,this.P,'EdgeAlpha',0)
%             title(fig_title)
            xlabel('X - axis')
            ylabel('Y - axis')
            zlabel('refracive index')
            colorbar
%             axis equal
        end
        
        function plot_gradient(this,figureNr)
            figure(figureNr)
            contour(this.X,this.Y,this.P)
            plot([0;0],[-1.2;1.2],'k')
            plot([-1.2;1.2],[0;0],'k')
            viscircles([0,0],1,'color','k')
            contour(this.X,this.Y,this.P)
            hold on
            %             quiver(this.X,this.Y,this.dX,this.dY)
            hold off
            axis equal
            %             title(fig_title)
            xlabel('X - axis')
            ylabel('Y - axis')
            zlabel('refracive index')
            %             colorbar
            axis equal
        end
        
        function plot_nIndex_line(this, figureNr)
            y = round(size(this.Y,1)/2);
            nline = this.P(y,y:end);
            x = this.X(y,y:end);
            figure(figureNr)
            hold on
            plot(x,nline)
        end
    end
end










