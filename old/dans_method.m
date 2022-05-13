function index_profile = dans_method(section, plot)
% n_points = section.n_shells;
index_profile = zeros(section.n_shells,2);

% Do initial refractive index aproximation.
for shell_ind = 2:section.n_shells
    if shell_ind == 2
        d = intersect_dist(section.phase_pos(shell_ind-1),section.phase_pos(shell_ind));
        index_profile(shell_ind,1) = section.n0 + section.phase_val(shell_ind)/d;
        continue
    end
    f = first_approx(shell_ind);
end

    function f = first_approx(shell_ind)
        
    end

    % function for finding distance through each shell
    function d = intersect_dist(r,x)
        if x >= r
            d = nan;
            return
        end
        d = 2*sqrt(r^2 - x^2);
    end

end