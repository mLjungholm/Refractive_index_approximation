clear
close all

% half_shift_pixels = [15 39 53 59 64 66];
shells = [66 64 59 53 39 15];   % Half shift in pixels
wavelength = 550*10^(-9);               % Unsertain value since it does not say in the paper.
radius = 13*10^(-6);                    % Radius for the crystaline cone
n0 = 1.34;                              % Refractive index of surounding medium
n_half = length(shells)-1;

phase_shift_pos = shells .* (radius/66);             % Phase shift positions in distance [m].
phase_shift_val = [0 (((1:1:n_half) .* wavelength) - wavelength/2)];


section = Section(phase_shift_val,phase_shift_pos);
section.linear_interp;
section.plot_phase_profile(1)


section.create_rays;
section.plot_shells(2)
section.plot_rays(2)
