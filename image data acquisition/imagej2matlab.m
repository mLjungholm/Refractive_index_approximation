% Function to convert from imageJ coordinates to matlab image coordinates.

% im_size in [x,y];
function p_list_out = imagej2matlab(p_list_in, im_size)
p_list_out = p_list_in;
p_nums = size(p_list_in,1);
for i = 1:p_nums
    p_list_out(i,2) = im_size(1)-p_list_in(i,2);
end
end