% Testscript running the old phase shift to refractive index algorithm. 


% Step one is to create the data needed for the script to run. This can
% either be done by importing an image and selecting lines as in the
% proposed app or simply a preset dataset.

% Manual imported data for the example in "Nilsson et.al 1983"

close all

%% Creating the phase shift profile

half_shift_pixels = [15 39 53 59 64 66];   % Half shift in pixels
wavelength = 550*10^(-9);               % Unsertain value since it does not say in the paper.
radius = 13*10^(-6);                    % Radius for the crystaline cone
n0 = 1.34;                              % Refractive index of surounding medium
n_half = length(half_shift_pixels);     % Number of measured halfshifts

phase_shift_pos = half_shift_pixels .* (radius/66);             % Phase shift positions in distance [m].
phase_shift_val = [(((n_half:-1:2) .* wavelength) - wavelength/2) 0]; % Phase shift values
% phase_shift_pos = phase_shift_pos.*10^6;
% phase_shift_val = phase_shift_val.*10^6;

figure(1)
hold on
plot(phase_shift_pos,phase_shift_val,'+')
xlabel('Radial position [m]')
ylabel('Phase shift value [m]')
title('Phase shift profile')

% Fitting the parabolic funtcion
% phase_shift_pos = [-fliplr(phase_shift_pos) phase_shift_pos];
% phase_shift_val = [fliplr(phase_shift_val) phase_shift_val];
f_profile = fit(phase_shift_pos',phase_shift_val','linearinterp');
% plot(f,phase_shift_pos(n_half+1:end),phase_shift_val(n_half+1:end))
plot(f_profile,phase_shift_pos,phase_shift_val)

%% Creating the ray list and refractive index matrix

% Divide the slice into shells.
n_shells = n_half-1;
shell_pos = zeros(n_shells,2);
shell_index = zeros(n_shells,2);


rays = Rays(n_shells,phase_shift_pos./10^6);

% ray_step = radius/10^3;
% 
% figure(2)
% title('Ray Trace')
% xlabel('[m]')
% hold on; axis equal
% viscircles([0,0],radius,'color','k','linewidth',0.5)
% for shell_ind = 1:n_shells
%     viscircles([0,0],rays_start(shell_ind,1),'color','k','linewidth',0.5)
% end
% plot([-radius; radius].*1.2, [0;0],'color','k')
% plot([0,0],[-radius; radius].*1.2,'color','k')
% 
% 
% function rays = set_up_rays(n_shells,start_positions)
% path = cell(n_shells,1);
% dir = zeros(n_shells,2);
% dir(:,1) = 0;
% dir(:,2) = -1;
% terminated = zeros(n_shells,1);
% pos = [start_positions', ones(n_shells,1).*max(start_positions)*1.3];
% rays.path = path;
% rays.dir = dir;
% rays.pos = pos;
% rays.terminated = terminated;
% rays.n_rays = n_shells;
% end

% function plot_rays(rays_pos,rays_dir)
% 
% end
