% v0 = [0 -1];
% p0 = [0.6 1.1];
% r0 = 1;
% ds = r0/10^3;
% 
% [p0,~,~] = circleIntersect(r0,p0,v0);
% 
% n0 = 1.3;
% n1 = 1.4;
% k = n0-n1;
% n = @(r) k*r + n1;
% 
% gV = -p0;
% r1 = norm(p0);
% ra = r1 + ds/2;
% rb = r1 - ds/2;
% na = n(ra);
% nb = n(rb);
% 
% a = dot(v0,gV);
% theta = acos(a);
% alpha = asin(na/nb*sin(theta));
% fi = alpha-theta;
% 
% r = ((na+nb)/2)/(sin(theta)*(nb-na)/abs(rb-ra));
% psi = -ds/r;
% % R = @(psi) [cos(psi) -sin(psi);
% %             sin(psi) cos(psi)];
% % v1 = R(psi)*v0';
% % v1 = v1';
% % or = p0' + r.*R(pi/2)*v0';
% % p1 = R(psi)*(p0'-or) + or;
% % p1 = p1';
% fprintf('meggit is :%f.6 rad \n',psi) 
% fprintf('snell in : %f.6 rad, snell out: %f.6 rad, diff is %f.6 rad \n',theta,alpha,fi)
% 
% 
% r = sqrt(x.^2 + y.^2);
% k = 0.6;
% dndr = @(x,y) k.*(k./sqrt(x.^2 + y.^2)).*[x y];
% dif = dndr(x',y');
% 
% ds = 0.01;
% xp = -1:ds:1;
% yp = -1:ds:1;
% [X,Y] = meshgrid(xp,yp);
% r =  sqrt(X.^2 + Y.^2);
% 
% n0 = 1.2;
% n1 = 2;
% k = n0-n1;
% n = k.*sqrt(X.^2 + Y.^2) + n1;
% % n(r > (1+ds/2) ) = n0;
% [dx,dy] = gradient(n,ds);
% 
% % figure(1)
% % hold on; axis equal; grid on
% % contour(X,Y,n)
% % quiver(X,Y,dx,dy)
% % plot([0;0],[-1.2;1.2],'k')
% % plot([-1.2;1.2],[0;0],'k')
% % plotCircle(1,2*pi,1)
% % xlabel('X - axis')
% % ylabel('Y - axis')
% 
% 
x = linspace(0,1/sqrt(2),100);
y = linspace(0,1/sqrt(2),100);
r = sqrt(x.^2 + y.^2);


n0 = 1.2;
n1 = 1.5;
k = n0-n1;
n = @(x,y) k.*sqrt(x.^2 + y.^2) + n1;
dndr = @(x,y) k./sqrt(x.^2 + y.^2).*[x y];
dif = dndr(x',y');
sqrt(dif(2,1)^2 + dif(2,2)^2)
% dif = dif.*n(x',y');
% 
% 
syms x y
ds = 0.1;
n = k.*sqrt(x^2 + y^2) + n1;
g = gradient(n,[x,y]);

[X, Y] = meshgrid(ds:ds:1,ds:ds:1);
G1 = subs(g(1),[x y],{X,Y});
G2 = subs(g(2),[x y],{X,Y});
quiver(X,Y,G1,G2)

