% class for holding information on the Sections. This includes number of
% shells, phase shift profile, refractive index and radius.

classdef Section < handle
    properties
%         radius
        phase_val
        phase_pos
        n_shells
        f_profile
        rays
        n0
%         shell_radius
    end
    
    methods
        function this = Section(phase_values,phase_positions)
            sorting_mat = [phase_values' phase_positions'];
            sorting_mat = sortrows(sorting_mat,2);
            sorting_mat = flipud(sorting_mat);
            this.phase_val = sorting_mat(:,1);
            this.phase_pos = sorting_mat(:,2);
            this.n_shells = length(phase_values);
        end
        
        function linear_interp(this)
            this.f_profile = fit(this.phase_pos,this.phase_val,'linearinterp');
        end
        
        function create_rays(this)
            this.rays = Rays(this.n_shells-1, this.phase_pos(2:end));
        end
        
        function plot_phase_profile(this,figure_nr)
            figure(figure_nr)
            hold on
            plot(this.phase_pos,this.phase_val,'+')
            xlabel('Radial position [m]')
            ylabel('Phase shift value [m]')
            title('Phase shift profile')
            plot(this.f_profile,this.phase_pos,this.phase_val)
        end
        
        function plot_shells(this,figure_nr)
            figure(figure_nr)
            title('Ray Trace')
            xlabel('[m]')
            hold on; axis equal
%             viscircles([0,0],radius,'color','k','linewidth',0.5)
            for shell_ind = 1:this.n_shells
                viscircles([0,0],this.phase_pos(shell_ind,1),'color','k','linewidth',0.5)
            end
            plot([-this.phase_pos(1); this.phase_pos(1)].*1.2, [0;0],'color','k')
            plot([0,0],[-this.phase_pos(1); this.phase_pos(1)].*1.2,'color','k')
        end
        
        function plot_rays(this, figure_nr)
            this.rays.plot_rays(figure_nr)
        end
    end
end