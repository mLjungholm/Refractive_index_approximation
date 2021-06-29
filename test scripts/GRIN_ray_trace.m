clear
close all

rayP = zeros(3,1);  % Ray point
rayV = zeros(3,1);  % Ray vector

dt = 0.1; % Stepsize 

[X,Y,P] = create_2d_grin(0.1,'parabolic','matrix');

[px,py] = gradient(P);
Dx = P.*px;
Dy = P.*py;

Dr = findDr(Dx,Dy,R,X,Y);
% A = dt



% Plot the contour lines and vectors in the same figure.
figure
contour(X,Y,P)
hold on
quiver(X,Y,Dx,Dy)
hold off
axis equal


function [rayP,rayV,flag] = findDr(rayP,rayV,Dx,Dy,X,Y)


end


function closestInd = findClosestInd(X,Y,rayP,rayV)