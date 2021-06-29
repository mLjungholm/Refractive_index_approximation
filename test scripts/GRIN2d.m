% circular 2D GRIN surface for testing GRIN ray trace
% The surface has a radius of 1 and is centered at [0,0]
classdef GRIN2d < handle
    
    properties
        X = []; % x coordinates
        Y = []; % y coordinates
        P = []; % index of refraction matrix
        gradient_type = ''; % Type of gradient
        data_format = '';   % data format, vector or matrix
        stepsize = [];      % grid spacing
        dX = [];            % nDx matrix. Dx = gradient of P(x)
        dY = [];            % nDy matrix. Dy = gradient of P(y)
        nEdge = 1;          % Refractive index att edge
        grinRange = [];
    end
    
    methods
        function this = GRIN2d(stepsize,gradient_type,data_format,refractiveIndex, grinRange)
            this.stepsize = stepsize;
            this.data_format = data_format;
            this.gradient_type = gradient_type;
            this.nEdge = refractiveIndex;
            this.grinRange = grinRange;
            xp = -1:stepsize:1;
            yp = -1:stepsize:1;
            [this.X,this.Y] = meshgrid(xp,yp);
            mask = sqrt(this.X.^2 + this.Y.^2);
            mask(mask>(1+stepsize/2)) = nan;
            mask(~isnan(mask)) = 1;
            
            switch gradient_type
                case 'linear'
                    % Linear space: n = 1.5 - r;  r[0,1]
                    nScale = this.grinRange(2)-this.grinRange(1);
                    this.P = (nScale.*sqrt(this.X.^2 + this.Y.^2)*-1 + this.grinRange(2)).*mask;
                    this.P(isnan(this.P)) = this.grinRange(1);
                case 'Circular'
                    % parabolic space: n = sqrt(1-r^2) + 1;  r[0,1]
                    nScale = this.grinRange(2)-this.grinRange(1);
                    r = sqrt(this.X.^2 + this.Y.^2).*mask;
                    this.P = sqrt(1 - r.^2).*nScale + this.grinRange(1);
                    this.P = real(this.P);
                    this.P(isnan(this.P)) = this.grinRange(1);
                case 'square'
                    % square space: n = 1.5;  r[0,1]
                    this.P = 1.5.*mask;
                    this.P(isnan(this.P)) = 1;
                case 'parabolic'
                    k = this.grinRange(1)-this.grinRange(2);
                    r = sqrt(this.X.^2 + this.Y.^2).*mask;
                    this.P = k.*r.^2 + this.grinRange(2);
                    this.P(isnan(this.P)) = this.grinRange(1);
                otherwise
                    disp('Error: Not valid type')
                    this.P = nan;
            end
            
            switch data_format
                case 'matrix'
                case 'vector'
                    this.X = reshape(this.X.*mask,[],1);
                    this.Y = reshape(this.Y.*mask,[],1);
                    this.P = reshape(this.P.*mask,[],1);
                    
                    this.X = this.X(~isnan(this.X));
                    this.Y = this.Y(~isnan(this.Y));
                    this.P = this.P(~isnan(this.P));
                otherwise
                    disp('Error: Not valid format')
                    this.P = nan;
                    this.X = nan;
                    this.Y = nan;
            end
            
            [px,py] = gradient(this.P,this.stepsize);
            this.dX = this.P.*px;
            this.dY = this.P.*py;
        end
        
        function plot_nIndex(this)
            fig_title = strcat('Refractive index -',this.gradient_type,' gradient');
            figure('Name','Refractive index plot','NumberTitle','off')
            surf(this.X,this.Y,this.P,'EdgeAlpha',0)
            title(fig_title)
            xlabel('X - axis')
            ylabel('Y - axis')
            zlabel('refracive index')
            colorbar
            axis equal
        end
        
        function plot_gradient(this,figureNr)
%             fig_title = strcat('Gradient -',this.gradient_type,' gradient');
%             figure('Name','Gradient plot','NumberTitle','off')
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










