% Class for holding information on each ray.
% Used in the old algorithm for calculating refractive index from phase
% shift profile.

classdef Rays < handle   
    properties
        n_rays;     % Number of rays in the bundle
        pos;        % Current position of each ray
        dir;        % Current direction of each ray
        path;       % Path for the rays. This will start as empty
        terminated; % =1 if ray has finished tracing.
        step_size;
    end
              
    methods
        function this = Rays(n_rays,shell_radius)
            this.n_rays = n_rays;
            this.pos = [shell_radius, ones(n_rays,1).*shell_radius(1)*1.2];
            this.path = cell(n_rays,1);
            this.dir = zeros(n_rays,2);
            this.dir(:,1) = 0;
            this.dir(:,2) = -1;
            this.terminated = zeros(n_rays,1);
            this.step_size = shell_radius(1)/10^3;
        end
        
        function plot_rays(this,figure_nr)
            figure(figure_nr)
            hold on
            if isempty(this.path{1})
                for ray_ind = 1:this.n_rays
                    path_x = [this.pos(ray_ind,1); (this.pos(ray_ind,1) + this.dir(ray_ind,1)*this.step_size*10^3*2)];
                    path_y = [this.pos(ray_ind,2); (this.pos(ray_ind,2) + this.dir(ray_ind,2)*this.step_size*10^3*2)];
                    plot(path_x,path_y,'color','r')
                end
            end
        end
    end
end