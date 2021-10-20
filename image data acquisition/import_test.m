% Trying to find a good way to import the data from the images

% ImageJ coordinates are ordered [x,y] with (0,0) in top left corner.
% Matlab image are ordered [y,x] with (0,0) in bottom left corner.
% I.e to transfer from imagej to matlab switch place on x,y and invert y.

close all

mask_im = uint8(zeros(size(test_cell_im)));
test_im = test_cell_im;

endy = size(test_cell_im,1);
% for i = 1:size(testarea,1)
%     mask_im(endy-testarea(i,2),testarea(i,1)) = 255;
%     test_im(endy-testarea(i,2),testarea(i,1)) = 255;
% end
for i = 1:size(mask_circle,1)
    mask_im(endy-mask_circle(i,2),mask_circle(i,1)) = 255;
    test_im(endy-mask_circle(i,2),mask_circle(i,1)) = 255;
end

pos = [2444,3293];
test_im = insertMarker(test_im,pos,'x','size',10);
imshow(test_im)