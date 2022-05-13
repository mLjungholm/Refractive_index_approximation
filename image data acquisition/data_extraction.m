
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this script for extracting the phase shift data from the images along
% the ROI lines.

% Each sett of data needs to include:
% - Start position
% - End position
% - Length [m] or pixel scale [m]
% - Pixel numbers
% - Pixel centers
% - Pixel values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initial data
end_points = zeros(2,2); % [x1, y1; x2, y2]
p_scale = 10^-6; % [m/pixel]
p_list = []; % Should contain all pixel coordinates in distance order (small -> large) change if not in order.
p_nums = size(pixel_list, 1); % Number of pixels
p_dist = zeros(pixel_nums,1); % Distance to center of each pixel from starting pixel.
p_val = zeros(p_nums,1); % List that contains the pixel value for each pixel.

% Extract data from image C
p_val(1) = C(p_list(1,2),p_list(1,1)); % Check the correct order of the indicies. XY coordinates in space is not the same as pixel coordinates.
% ImageJ coordinates are ordered [x,y] with (0,0) in top left corner.
% Matlab image are ordered [y,x] with (0,0) in bottom left corner.
% I.e to transfer from imagej to matlab switch place on x,y and invert y.

for pixel_ind = 2:pixels_nums  % This could be made in a parallel manner. Nut shure if it is nessecary though. Maybe fast enough.
    
end