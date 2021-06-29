function [X,Y,P] = create_2d_grin(stepsize, type, format)

% Creating grid
% stepsize = 0.01;
xp = -1:stepsize:1;
yp = -1:stepsize:1;
[X,Y] = meshgrid(xp,yp);
mask = sqrt(X.^2 + Y.^2);
mask(mask>1) = nan;
mask(~isnan(mask)) = 1;
% x = reshape(X.*mask,[],1);
% y = reshape(Y.*mask,[],1);

switch type

    % Linear space: n = 1.5 - r;  r[0,1]
    case 'linear'
        P = (1.5 - sqrt(X.^2 + Y.^2)./2).*mask;

    % parabolic space: n = sqrt(1-r^2) + 1;  r[0,1]
    case 'circular'
        P = sqrt(1-sqrt(X.^2 + Y.^2).^2.*mask)./2+1;
        
    % square space: n = 1.5;  r[0,1]
    case 'square'
        P = 1.5.*mask;
    case 'parabolic'
        P = sqrt(X.^2 + Y.^2).^2;
    otherwise
        disp('Error: Not valid type')
        P = nan;
end

switch format
    case 'matrix'
    case 'vector'
        X = reshape(X.*mask,[],1);
        Y = reshape(Y.*mask,[],1);
        P = reshape(P.*mask,[],1);
        
        X = X(~isnan(X));
        Y = Y(~isnan(Y));
        P = P(~isnan(P));
    otherwise
        disp('Error: Not valid format')
        P = nan;
        X = nan;
        Y = nan;
end

end


%% Visualization

% figure(2)
% nl(isnan(nl))=1;
% surf(X,Y,nl,'EdgeAlpha',0)
% title('Refractive index gradient  -  Linear [1.5:1]')
% xlabel('X - axis')
% ylabel('Y - axis')
% zlabel('refracive index')
% colorbar
% axis equal
% 
% 
% figure(3)
% surf(X,Y,np,'EdgeAlpha',0)
% np(isnan(np))=1;
% title('Refractive index gradient  -  parabolic [1.5:1]')
% xlabel('X - axis')
% ylabel('Y - axis')
% zlabel('refracive index')
% colorbar
% axis equal
% 
% 
% figure(4)
% ns(isnan(ns))=1;
% surf(X,Y,ns,'EdgeAlpha',0)
% title('Refractive index gradient  -  tophat [1.5]')
% xlabel('X - axis')
% ylabel('Y - axis')
% zlabel('refracive index')
% colorbar
% axis equal
