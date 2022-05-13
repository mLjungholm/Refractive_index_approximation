% Script version of the graded index refraction aproximation algorithm

% v.0.1

% By: Mikael Ljungholm 2020

% The script takes interferometric images as inputs and calculates the
% graded refractive index of a rotationaly symetric object.

%% Setp 1. Load image
% 
% [filename, path] = uigetfile('*');
% filePath = strcat(path,filename);
% I = imread(filePath);


%% Step 2. Create ROIs

figure(1)
imshow(I)
roi_list = myLoopingFcn;

disp('Draw Lines by click and hold left mouse button')
disp('Press any key to cansel loop')


function roi_list = myLoopingFcn() 
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf
set(gcf, 'KeyPressFcn', @myKeyPressFcn)
n_roi = 0;
roi_list = {};
while ~KEY_IS_PRESSED
    drawnow
    n_roi = n_roi + 1;
    roi_list{n_roi} = drawline('InteractionsAllowed','none');
end
% disp('loop ended')
end

function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
disp('Draw final line')
end