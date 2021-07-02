% clear
% close all

% Create GRIN
grinRange = [1.3,1.45];
grin = GRIN2d(0.001,'parabolic','matrix',1.3,grinRange);
% grin.plot_nIndex
% % Create source
start = [0,2];
v = [0,-1];
rayNr = 20;
sourceWidth = 2.2;
n_source = 1.3;
s2 = source2d(start,v, rayNr, sourceWidth, n_source,'half');

% % Testplot source and grin
% figure(2)
% hold on, axis equal, grid on
% plotCircle(1,2)
% grin.plot_gradient(1)
% s.plotP(1)
tic
% halfstop = true;
stopLine = -1;
rayTrace2dGRIN_parallel(s2,grin,0.001,stopLine);
% rayTrace2d_homogenious(s)
% % rayTrace2dGRIN(s,grin,0.001,halfstop);
toc
s2.plotTrace(1,'b')
% s.plotBacktrack(0,2)
% plot([0,0],[1.2,-1.2],'k')
% totalPath = s.getTotalPath;
% s.pathHist(2)


%% This is afunctioning test ray trace for a single ray
% % GRIN
% [X,Y,P] = create_2d_grin(0.01, 'linear', 'matrix');
% [px,py] = gradient(P,0.01);
% dX = P.*px;
% dY = P.*py;
% 
% % Ray
% V = [1; 0];
% R = [-2; 0.5];
% n1 = 1;
% n2 = 1;
% 
% [p,intersect] = circleLineIntersect(1,R',V');
% 
% if ~intersect
%     return
% end
% 
% dt = 0.01;
% ended = 0;
% maxSteps = round(2/dt*1.5);
% stepNr = 3;
% rayPath = zeros(maxSteps,2);
% rayPath(:,1) = nan; rayPath(:,2) = nan;
% rayVpath = rayPath;
% raydV = rayPath;
% rayPath(1,:) = R';
% rayPath(2,:) = p;
% R = p';
% rayVpath(1,:) = V';
% N = -p./sqrt(p(1)^2+p(2)^2);
% [vNew, reflected] = Snell(V', N, n1, n2);
% rayVpath(2,:) = vNew;
% V = vNew';
% inside = 0;
% while ~inside
%     closest = findClosest(X,Y,R);
%     dP = [dX(closest(2),closest(1));dY(closest(2),closest(1))];
%     if isnan(dP(1)) || isnan(dP(2))
%         R = R + dt.*V;
%         rayPath(stepNr,:) = R';
%         rayVpath(stepNr,:) = V;
%         stepNr = stepNr + 1;
%     else
%         inside = 1;
%     end
% end
%     
% while ~ended
%     closest = findClosest(X,Y,R);
%     dP = [dX(closest(2),closest(1));dY(closest(2),closest(1))];
%     raydV(stepNr,:) = dP';
%     A = dt*dP;
%     Bv = R + (dt/2)*V + 1/8*dt*A;
%     closest = findClosest(X,Y,Bv);
%     B = dt*[dX(closest(2),closest(1));dY(closest(2),closest(1))];
%     Cv = R + dt*V + 1/2*dt*B;
%     closest = findClosest(X,Y,Cv);
%     C = dt*[dX(closest(2),closest(1));dY(closest(2),closest(1))];
%     R2 = R + dt*(V + 1/6*(A + 2*B));
%     V2 = V + 1/6*(A + 4*B + C);
%     if isnan(R2(1)) || isnan(R2(2))
%         ended = 1;
%         break
%     end   
%     rayPath(stepNr,:) = R2';
%     rayVpath(stepNr,:) = V';
%     R = R2;
%     V = V2;
%     stepNr = stepNr + 1;
%     if stepNr >= maxSteps
%         ended = 1;
%         break
%     end
% end
% 
% % figure(1)
% % xlabel('X - axis')
% % ylabel('Y - axis')
% % axis equal
% % hold on
% % contour(X,Y,P)
% % plot(R(1),R(2),'.')
% % quiver(R(1),R(2),V(1),V(2),0.5,'b')
% % plot([0;0],[-1.2;1.2],'k')
% % plot([-1.2;1.2],[0;0],'k')
% % viscircles([0,0],1,'color','k')
% % quiver(R(1),R(2),dP(1),dP(2),50,'r')
% % grid on
% 
% figure(1)
% xlabel('X - axis')
% ylabel('Y - axis')
% axis equal
% hold on
% contour(X,Y,P)
% plot([0;0],[-1.2;1.2],'k')
% plot([-1.2;1.2],[0;0],'k')
% viscircles([0,0],1,'color','k')
% plot(rayPath(:,1),rayPath(:,2),'r')
% % quiver(rayPath(:,1),rayPath(:,2),raydV(:,1),raydV(:,2),1,'r')
% grid on
% 
% 
% function closest = findClosest(X,Y,p)
% [~,closestX] = min(abs(X(1,:)-p(1)));
% [~,closestY] = min(abs(Y(:,1)-p(2)));
% closest = [closestX,closestY];
% end
% 
% function [p,intersect] = circleLineIntersect(r,p0,v)
%  L = -p0;
%  tca = v(1)*L(1) + v(2)*L(2);
%  if tca < 0
%      intersect = 0;
%      p = nan;
%      return
%  end
%  d = sqrt(L(1)^2 + L(2)^2 - tca^2);
%  if d < 0 || d > r
%      intersect = 0;
%      p = nan;
%      return
%  end
%  thc = sqrt(r^2-d^2);
%  t0 = tca - thc;
%  p = p0 + t0*v;
%  intersect = 1;
% end











