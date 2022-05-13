% Improved and simplified version of Dan's method for calculating
% refractive-index gradient from phase-shift profiles

close all
clear

% Measured values
shells = [66 64 59 53 39 15];       % Half shift in pixels
wavelength = 550*10^(-9);           % Unsertain value since it does not say in the paper.
radius = 13*10^(-6);                % Radius for the crystaline cone
n0 = 1.34;                          % Refractive index of surounding medium

shell_r = shells .* (radius/66);    % Phase shift positions in distance [m].
phase_shift_val = [0 (((1:1:n_half) .* wavelength) - wavelength/2)];








