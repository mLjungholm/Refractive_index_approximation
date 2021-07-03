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

n0 = 1.3;
n1 = 1.4;
rn = 1;

k = (n1-n0)/rn;
n = @(x,y) k.*sqrt(rn^2 - x.^2 - y.^2) + n0;
dndr = @(r) -r*k/sqrt(rn^2 - r^2);

x = linspace(0,1/sqrt(2),100)';
y = linspace(0,1/sqrt(2),100)';

plot(sqrt(2.*x.^2),n(x,y))
        
        