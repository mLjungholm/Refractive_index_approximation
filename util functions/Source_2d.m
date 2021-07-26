% In use

classdef Source_2d < handle
    % Parent class for the source objects
    properties
        source_type = '2D Source';
        def = 'Undefined';
        id = 'Undefined';
        % Souce variables:
        nRays = 0;  % Number of rays in source
        terminated = []; % Termination vector
        nStart = 1; % Refractive index of start medium
        
        % Ray variables:
        P = [];             % Current point vector
        V = [];             % Current direction vector
        p0 = []; v0 = [];   % Start vector and coordinates for each ray
        path = {};          % Path vector
        totalPath = [];     % Total paht traveled
        phase = [];         % Acumilated phase shift.
        projection = [];    % Ray projection ontwo a line
        projection_base_v   % Base vector for the projections
        projection_center   % Centerpoint for the projection
        
        % Test variables:
        diviation = {}; % Vector diviation in each step [rad]
        stepError = {}; % Step length differnce from 'ds' in each step
    end
    
    methods
        function this = Source_2d(startP,startV,varargin)
            size_point_in = size(startP);
            if ~isempty(varargin)
                %             if size_point_in(1) == 1 || size_point_in(2) == 1
                try validateattributes(varargin{1},{'numeric'},{'scalar','real','integer'})
                catch
                    fprintf('Error: Number of rays = varargin{1} needs to be a integer value \n')
                    return
                end
                try validateattributes(varargin{2},{'numeric'},{'scalar','real'})
                catch
                    fprintf('Error: Source width = varargin{2} needs to be a numeric value \n')
                    return
                end
                this.def = 'Equal ray spacing';
                this.nRays = varargin{1};
                %                 this.sWidth = varargin{2};
                this.terminated = zeros(this.nRays,1);
                this.path = cell(this.nRays,1);
                this.phase = zeros(this.nRays,1);
                this.totalPath = zeros(this.nRays,1);
                this.diviation = cell(this.nRays,1);
                this.stepError = cell(this.nRays,1);
                this.V = ones(this.nRays,2).*startV;
                this.v0 = this.V;
                this.P = ones(this.nRays,2).*startP;
                % Rotate source line to match starting vector
                if this.nRays > 1
                    rg = linspace(-varargin{2},varargin{2},this.nRays);
                else
                    rg = 0;
                end
                this.P(:,2) = this.P(:,2)+rg';
                theta = -acos(startV(1)./sqrt(startV(1)^2+startV(2)^2));
                this.projection = zeros(this.nRays,2);
                if startV(2) < 0
                    theta = -theta;
                end
                if theta ~= 0
                    this.P = this.P - startP;
                    rot_theta = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                    this.P = this.P*rot_theta;
                    this.P = this.P + startP;
                    this.p0 = this.P;
                end
            else
                this.def = 'Defined ray positions';
                this.nRays = size(startP,1);
                this.P = startP;
                this.p0 = this.P;
                this.V = ones(this.nRays,2).*startV;
                this.v0 = this.V;
                this.terminated = zeros(this.nRays,1);
                this.path = cell(this.nRays,1);
                this.phase = zeros(this.nRays,1);
                this.totalPath = zeros(this.nRays,1);
                this.projection = zeros(this.nRays,2);
                this.diviation = cell(this.nRays,1);
                this.stepError = cell(this.nRays,1);
            end
        end
        
        % Resets the source to starting perameters
        function resetSource(this)
            this.P = this.v0;
            this.V = this.v0;
            this.terminated = zeros(this.nRays,1);
            this.path = cell(this.nRays,1);
            this.phase = zeros(this.nRays,1);
            this.totalPath = zeros(this.nRays,1);
            this.diviation = cell(this.nRays,1);
            this.stepError = cell(this.nRays,1);
            this.projection = zeros(this.nRays,2);
        end
        
        % Plotting rays
        function plotRays(this,color,varargin)
            if isempty(varargin)
                if isempty(this.nRays)
                    fprintf('Error in source.plotRays (%s)\n',this.id)
                    fprintf('Rays are not defined \n\n')
                    return
                end
                hold on
                if  any(cellfun(@isempty,this.path))
                    quiver(this.P(:,1),this.P(:,2),this.V(:,1),this.V(:,2),color)
                else
                    for rayNr = 1:this.nRays
                        plot(this.path{rayNr}(:,1),this.path{rayNr}(:,2),color)
                    end
                end
            else
                try
                    validateattributes(varargin{1},{'numeric'},{'2d','real','integer'})
                    rayInds = varargin{1};
                    if isempty(this.nRays)
                        fprintf('Error in source.plotRays (%s)\n',this.id)
                        fprintf('Rays are not defined \n\n')
                        return
                    end
                    hold on
                    for rayInd = 1:length(rayInds)
                        id = rayInds(rayInd);
                        if isempty(this.path{1})
                            quiver(this.P(id),this.P(id,2),this.V(id,1),this.V(id,2),color)
                        else
                            plot(this.path{id}(:,1),this.path{id}(:,2),color)
                        end
                    end
                catch
                    fprintf('Error in source.plotRays (%s)\n',this.id)
                    fprintf('Input variable needs to be numerical indices of the rays or empty \n\n')
                end
            end
        end
        
        % Project the rays onto line defined by a point "p" and vector "v".
        % Returns "nan" if there is no intersection and "inf" if the ray
        % lies on the defined line
        function intersection = projectRays(this,p,v,varargin)
            this.projection_base_v = v;
            this.projection_center = p;
            intersection = zeros(this.nRays,2);
            if ~isempty(varargin)
                if isequal(varargin{1},'back') || varargin{1} == -1
                    c = 1;
                end
            else
                c = -1;
            end
            for rayInd = 1:this.nRays
                A = [this.V(rayInd,:)' c.*v'];
                x = p' - this.P(rayInd,:)';
                is = A\x;
                if isnan(is(1))
                    intersection(rayInd,:) = [nan nan];
                    %                 elseif is(1) < 1
                    %                     intersection(rayInd,:) = [nan nan];
                else
                    pi = p' + is(2).*v';
                    %                     px = p(1) + is(1).*v(1);
                    %                     py = p(2) + is(2).*v(2);
                    %                     pi = [px ,py];
                    vi = -c.*(pi' - this.P(rayInd,:));
                    vi = vi./norm(vi);
                    dir = dot(this.V(rayInd,:),vi);
                    if dir < 0
                        intersection(rayInd,:) = [nan nan];
                    else
                        intersection(rayInd,:) = -pi';
                    end
                end
            end
            this.projection = intersection;
        end
        
        % Plot projected rays
        function plotProjection(this,color,varargin)
            if isempty(this.projection)
                fprintf('Error in source.plotProjection (%s)\n',this.id')
                fprintf('No projection defined \n\n')
                return
            elseif isempty(varargin)
                if isempty(this.nRays)
                    
                    fprintf('Error in source.plotProjection (%s)\n',this.id')
                    fprintf('Rays are not defined \n\n')
                    return
                end
                hold on
                for rayNr = 1:this.nRays
                    px = [this.P(rayNr,1);this.projection(rayNr,1)];
                    py = [this.P(rayNr,2);this.projection(rayNr,2)];
                    plot(px,py,color)
                end
            else
                try
                    validateattributes(varargin{1},{'numeric'},{'2d','real','integer'})
                    rayInds = varargin{1};
                    if isempty(this.nRays)
                        
                        fprintf('Error in source.plotProjection (%s)\n',this.id')
                        fprintf('Rays are not defined \n\n')
                        return
                    end
                    hold on
                    for id = 1:length(rayInds)
                        plotLine(this.P(rayInds(id),:),this.projection(rayInds(id),:),color)
                    end
                catch
                    
                    fprintf('Error in source.plotProjection (%s)\n',this.id')
                    fprintf('Input variable needs to be numerical indices of the rays or empty \n\n')
                end
            end
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
