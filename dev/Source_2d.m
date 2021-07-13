classdef Source_2d < handle
    % Parent class for the source objects
    properties
        source_type = '2D Source';
        def = 'Undefined';
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
        
        % Test variables:
        diviation = {}; % Vector diviation in each step [rad]
        stepError = {}; % Step length differnce from 'ds' in each step
    end
    
    methods
        function this = Source_2d(startP,startV,varargin)
            size_point_in = size(startP);    
            if size_point_in(1) == 1 || size_point_in(2) == 1
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
                    rg = linspace(-varargin{2},varargin{2},nRays);
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
        function plotRays(this,varargin)
            if nargin == 0
                if isempty(this.nRays)
                    fprintf('Error: Rays are not defined \n')
                    return
                end
                hold on
                if isempty(this.path{1})
                    quiver(this.P(:,1),this.P(:,2),this.V(:,1),this.V(:,2),'b')
                else
                    for rayNr = 1:this.nRays
                        plot(this.path{rayNr}(:,1),this.path{rayNr}(:,1),'r')
                    end
                end
            else
                try 
                    validateattributes(varargin{1},{'numeric'},{'2d','real','integer'})
                    rayInds = varargin{1};
                    if isempty(this.nRays)
                        fprintf('Error: Rays are not defined \n')
                    return
                    end
                    hold on
                    for rayInd = 1:length(rayInds)
                        id = rayInds(rayInd);
                        if isempty(this.path{1})
                            quiver(this.P(id),this.P(id,2),this.V(id,1),this.V(id,2),'b')
                        else
                            plot(this.path{id}(:,1),this.path{id}(:,1),'r')
                        end                        
                    end
                catch
                    fprintf('Error: Input variable needs to be numerical indices of the rays or empty /n')
                end
            end
        end
        
        % Project the rays onto line defined by a point "p" and vector "v".
        % Returns "nan" if there is no intersection and "inf" if the ray
        % lies on the defined line
        function intersection = projectRays(this,p,v,varargin)
            intersection = zeros(this.nRays,2);
            if isequal(varargin{1},'back') || varargin{1} == -1
                c = 1;
            else
                c = -1;
            end
            for rayInd = 1:this.nRays
                A = [this.V(rayInd,:)' c.*v'];
                x = p' - this.P(rayInd,:)';
                is = A\x;
                if isnan(is(1))
                    intersection(rayInd,:) = [nan nan];
                elseif is(1) < 1
                    intersection(rayInd,:) = [nan nan];
                else
                    intersection(rayInd,:) = p + is(2).*v;
                end
            end
            this.projection = intersection;
        end
        
        % Plot projected rays
        function plotProjection(this,varargin)
            if isempty(this.projection)
                    fprintf('Error: No projection defined \n')
                    return
            elseif nargin == 0
                if isempty(this.nRays)
                    fprintf('Error: Rays are not defined \n')
                    return
                end
                hold on
                for rayNr = 1:this.nRays
                    plotLine(this.P(rayNr,:),this.projection(rayNr,:),'--r')
                end
            else
                try 
                    validateattributes(varargin{1},{'numeric'},{'2d','real','integer'})
                    rayInds = varargin{1};
                    if isempty(this.nRays)
                        fprintf('Error: Rays are not defined \n')
                    return
                    end
                    hold on
                    for id = 1:length(rayInds)
                        plotLine(this.P(rayInds(id),:),this.projection(rayInds(id),:),'--r')
                    end
                catch
                    fprintf('Error: Input variable needs to be numerical indices of the rays or empty /n')
                end
            end
        end
    end
end
    