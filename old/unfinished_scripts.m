        function add_phase_values(this,phase_vals,phase_pos)
            if size(phase_vals,1) < size(phase_vals,2)
                this.T = [this.T; phase_vals'];
            else
                this.T = [this.T; phase_vals];
            end
            if size(phase_pos,1) < size(phase_pos,2)
                this.shell_r = [this.shell_r; phase_pos'];
            else
                this.shell_r = [this.shell_r; phase_pos];
            end
            sortT = sortrows([this.T this.shell_r],2);
            this.T = sortT(:,1);
            this.shell_r = sortT(:,2);
            this.rays_x = this.shell_r(2:end);
            this.n_rays = 2
        end