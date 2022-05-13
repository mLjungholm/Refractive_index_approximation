% Trying to find a good way to import the data from the images

% ImageJ coordinates are ordered [x,y] with (0,0) in top left corner.
% Matlab image are ordered [y,x] with (0,0) in bottom left corner.
% I.e to transfer from imagej to matlab switch place on x,y and invert y.

% close all

% mask_im = uint8(zeros(size(test_cell_im)));
% test_im = test_cell_im;
% mask2 = imagej2matlab(mask_circle,size(test_im));
% % mask2 = mask_circle;
% 
% endy = size(test_cell_im,1);
% 
% for i = 1:size(mask_circle,1)
%     test_im(mask2(i,2),mask2(i,1)) = 255;
% end
% 
% pos = [2444,3293];
% test_im = insertMarker(test_im,pos,'x','size',10);
% figure(1)
% imshow(test_im)

%% Try to fit a spline on the data points


pos = [549.65295,683.24347];
m = [28.865107,199.44662];
% m2 = 

p = pos + m;