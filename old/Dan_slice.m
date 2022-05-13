classdef Dan_slice < handle
    properties
        T = [];         % The measured phase shift values
        shell_r = [];   % Radius for the shells and also the x-coordinates for the phase-shift values.
        rays_x = [];    % x coordinate for the rays ray_x(i) = shell(i+1).
        n0 = 0;         % Refractive index of the surounding medium.
        n = [];         % Refractive index for each shell.
        n_rays = 0;     % Number of rays.
        raypaths = [];  % Path coordinates for each ray.
    end
    
    methods
        %Instantiation function. The phase-shift values and positions
        %should be in the form of culomn vectors.
        function this = Dan_slice(phase_shift_vals,phase_shift_pos,n0)
            this.T = phase_shift_vals;
            this.shell_r = phase_shift_pos;
            this.n0 = n0;
            this.rays_x = phase_shift_pos(2:end);
            this.n_rays = length(this.rays_x);
            this.n = zeros(this.n_rays,1);
            this.raypaths = cell(this.n_rays,1);
        end
        
        function approximate_n(this, step_size)
            this.first_approx;
            for ray_ind = 1:this.n_rays
                trace_ray(ray_ind, step_size);
            end
        end
        
        % Function for doing the first approximation of each shell
        % refractive-index.
        function first_approx(this)
            this.n(1) = this.n0 + this.T(1)/(2*sqrt(this.shell_r(1)^2 - this.rays_x(1)^2));
            for ray_i = 2:this.n_rays
                Tt = zeros(ray_i-1,1);
                dt = zeros(ray_i,1);
                for shell_i = 1:ray_i
                    dt(shell_i) = 2*sqrt(this.shell_r(shell_i)^2 - this.rays_x(ray_i));
                end
                for shell_i = (length(shell_i)-1):-1:1
                    dt(shell_i) = dt(shell_i) - dt(shell_i + 1);
                end
                for shell_i = 1:ray_i-1
                    Tt(shell_i) = dt(shell_i)*(this.n(shell_i) - this.n0);
                end
                this.n{ray_i+1} = (this.T(ray_i) - sum(Tt))/dt(end) + this.n0;
            end
        end
        
        % Main ray-tracing algorithm. 
        function trace_ray(this, ray_i, step_size)
            % Setup raypath and starting position
            raypath = zeros(round(this.shell_r(1)/step_size),2);
            pos = [this.ray_x(ray_i) sqrt(this.shell_r(1)^2-this.ray_x(ray_i)^2)];
            raypath(1,:) = [this.ray_x(ray_i) this.shell_r(1)*1.2];
            raypath(2,:) = pos;
                    
            % Trace ray and calulate refractive index untill threshold is
            % met.
            acc_threshold = 0;
            while ~acc_threshold
                
                no = this.n0;
                nc = this.n(1);
                current_shell = 2;
                % Trace ray untill it emerges from the shells
                emerged = 0;
                while ~emerged
                    % trace one step.
                    D = this.shell_r(current_shell-1) - this.shell_r(current_shell); % Distance between the current shells
                    vn = point./sqrt(point(1)^2 + point(2)^2);  % Vecor normal of the refractive-index gradient
                    theta = acosd(dir(1)*vn(1) + dir(2)*vn(2)); % 
                    r = ((nc+no)/2) / (sind(theta)*(no-nc) / D); % Local curvature or radius
                    
                    % Check if ray has exit system
                    if sqrt(pos(1)^2 + pos(2)^2) >= this.shell_r(1)
                        emerged = 1;
                    end
                end
                % Calculate new refractive index
            end            
            this.raypaths{ray_i} = raypath;
        end
        
        
        function plot_shells(this,figure_nr)
            figure(figure_nr)
            title('Ray Trace')
            xlabel('[m]')
            hold on; axis equal
            for shell_ind = 1:length(this.shell_r)
                viscircles([0,0],this.shell_r(shell_ind),'color','k','linewidth',0.5)
            end
            plot([-this.shell_r(1); this.shell_r(1)].*1.2, [0;0],'color','k')
            plot([0,0],[-this.shell_r(1); this.shell_r(1)].*1.2,'color','k')
        end
        
        function plot_rays(this, figure_nr)
            figure(figure_nr)
            hold on
            for ray_i = 1:this.n_rays
                path_x = [this.rays_x(ray_ind,1); this.rays_x(ray_ind,1)];
                path_y = [this.shell_r(1); -this.shell_r(1)];
                plot(path_x,path_y,'color','r')
            end
        end
    end
end