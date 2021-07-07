n0 = 1;
n1 = 1.5;
k = n0-n1;
dt = 0.01;
v0 = [0 -1];

x = -1:ds:1;
y = 1:-ds:-1;
inds = length(x);

psi_rk = zeros(inds);

for xi = 1:inds
    for yi = 1:inds
        psi_rk(yi,xi) = getDeflection_RK(x(xi),y(yi),v0,ds,n0,n1);
    end
end

close all
figure(1)
imagesc(psi_rk)
title('Deflection [rads] Runge-Kutta')

ds = 0.1;
x = -1:ds:1;
y = 1:-ds:-1;
inds = length(x);

% psi_mg = getDeflection_MG(1,0,v0,ds,n0,n1);

psi_mg = zeros(inds);
for xi = 1:inds
    for yi = 1:inds
        psi_mg(yi,xi) = getDeflection_MG(x(xi),y(yi),v0,ds,n0,n1);
    end
end

figure(2)
imagesc(psi_mg)
title('Deflection [rads] Meggit')



function psi = getDeflection_RK(x,y,v0,dt,n0,n1)
k = n0-n1;
n = @(x,y) k.*sqrt(x.^2 + y.^2) + n1;
nDivn = @(R) (k.*sqrt(R(1).^2 + R(2).^2) + n1).*k./sqrt(R(1).^2 + R(2).^2).*[R(1) R(2)];
T = v0.*n(x,y);
A1 = dt.*nDivn([x,y]);
A2 = dt.*nDivn([x,y] + dt/2.*T +dt/8.*A1);
A3 = dt.*nDivn([x,y] + dt.*T +dt/2.*A2);
T1 = T + 1/6.*(A1 + 4.*A2 + A3);
v1 = T1./n(x,y);
psi = acos((v0(1)*v1(1) + v0(2)*v1(2))./norm(v1)./norm(v0));
end


function psi = getDeflection_MG(x,y,v0,ds,n0,n1)
k = n0-n1;
r = sqrt(x.^2 + y.^2);
n = @(x,y) k.*r + n1;
dndr = k;
gv = -[x,y];
theta = acos((gv(1)*v0(1) + gv(2)*v0(2))./(norm(gv)*norm(v0)));
psi = ds*sin(theta)*dndr/n(r);
psi = abs(psi);
end
% close all
% figure(1)
% axis equal;
% % quiver(X,Y,Dx,Dy)
% imagesc(absD)
% title('abs(n*grad(n)) RK')
% 
% figure(2)
% axis equal;
% imagesc(absDr)
% title('abs(grad(n)) RK') 
% 
% figure(3)
% axis equal;
% quiver(X,Y,dndx,dndy)
% title('grad(n) RK')
% 
% V = [0 -1];
% D = @(x,y) (k./sqrt(x.^2 + y.^2)).*(k.*sqrt(x.^2 + y.^2) + n1);
% A1x = ds.*D(X,Y).*X;
% A1y = ds.*D(X,Y).*Y;
% 
% X2 = X + ds/2.*V(1) + ds/8.*A1x;
% Y2 = Y + ds/2.*V(2) + ds/8.*A1y;
% A2x = ds.*D(X2,Y2).*X2;
% A2y = ds.*D(X2,Y2).*Y2;
% 
% X3 = X + ds.*V(1) + ds/2.*A2x;
% Y3 = Y + ds.*V(2) + ds/2.*A2y;
% A3x = ds.*D(X3,Y3).*X3;
% A3y = ds.*D(X3,Y3).*Y3;
% 
% Vkx = V(1) + 1/6.*(A1x + 4.*A2x + A3x);
% Vky = V(2) + 1/6.*(A1y + 4.*A2y + A3y);
% absV = sqrt(Vkx.^2 + Vky.^2);
% psi = acos((Vkx.*V(1) + Vky.*V(2))./absV);
% 
% figure(4)
% axis equal;
% % quiver(X,Y,Vkx,Vky)
% imagesc(psi)
% title('deflection RK')
% 
% 
% ds=0.1;
% [X, Y] = meshgrid(-1:ds:1,-1:ds:1);
% r = sqrt(X.^2 + Y.^2);
% n = k.*sqrt(X.^2 + Y.^2) + n1;
% 
% dndr = X./X.*k;
% dndr(isnan(dndr)) = k; 
% Rx = -X;
% Ry = -Y;
% theta = acos((Rx.*V(1) + Ry.*V(2))./r);
% psi2 = ds.*sin(theta).*dndr./n;
% 
% figure(5)
% axis equal;
% % quiver(X,Y,Vkx,Vky)
% imagesc(psi2)
% title('deflection Meggit')
% 
% 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% n0 = 1;
% n1 = 1.5;
% k = n0-n1;
% % n = @(x,y) k.*sqrt(x.^2 + y.^2) + n1;
% % dndr = @(x,y) k./sqrt(x.^2 + y.^2).*[x y];
% % dif = dndr(x',y');
% % sqrt(dif(2,1)^2 + dif(2,2)^2)
% % dif = dif.*n(x',y');
% 
% 
% 
% %linear  n(r)= k*r + n1
% 
% ds=0.1;
% [X, Y] = meshgrid(-1:ds:1,-1:ds:1);
% r = sqrt(X.^2 + Y.^2);
% n = k.*sqrt(X.^2 + Y.^2) + n1;
% dndr = k./sqrt(X.^2 + Y.^2);
% dndy = dndr.*Y;
% dndx = dndr.*X;
% absDr = sqrt(dndx.^2 + dndy.^2); absDr(isnan(absDr)) = max(max(absDr)); 
% D = k.*(k + n1./sqrt(X.^2 + Y.^2));
% Dy = D.*Y;
% Dx = D.*X;
% absD = sqrt(Dx.^2 + Dy.^2); absD(isnan(absD)) = max(max(absD)); 
% 
% close all
% figure(1)
% axis equal;
% % quiver(X,Y,Dx,Dy)
% imagesc(absD)
% title('abs(n*grad(n)) RK')
% 
% figure(2)
% axis equal;
% imagesc(absDr)
% title('abs(grad(n)) RK') 
% 
% figure(3)
% axis equal;
% quiver(X,Y,dndx,dndy)
% title('grad(n) RK')
% 
% V = [0 -1];
% D = @(x,y) (k./sqrt(x.^2 + y.^2)).*(k.*sqrt(x.^2 + y.^2) + n1);
% A1x = ds.*D(X,Y).*X;
% A1y = ds.*D(X,Y).*Y;
% 
% X2 = X + ds/2.*V(1) + ds/8.*A1x;
% Y2 = Y + ds/2.*V(2) + ds/8.*A1y;
% A2x = ds.*D(X2,Y2).*X2;
% A2y = ds.*D(X2,Y2).*Y2;
% 
% X3 = X + ds.*V(1) + ds/2.*A2x;
% Y3 = Y + ds.*V(2) + ds/2.*A2y;
% A3x = ds.*D(X3,Y3).*X3;
% A3y = ds.*D(X3,Y3).*Y3;
% 
% Vkx = V(1) + 1/6.*(A1x + 4.*A2x + A3x);
% Vky = V(2) + 1/6.*(A1y + 4.*A2y + A3y);
% absV = sqrt(Vkx.^2 + Vky.^2);
% psi = acos((Vkx.*V(1) + Vky.*V(2))./absV);
% 
% figure(4)
% axis equal;
% % quiver(X,Y,Vkx,Vky)
% imagesc(psi)
% title('deflection RK')
% 
% 
% ds=0.1;
% [X, Y] = meshgrid(-1:ds:1,-1:ds:1);
% r = sqrt(X.^2 + Y.^2);
% n = k.*sqrt(X.^2 + Y.^2) + n1;
% 
% dndr = X./X.*k;
% dndr(isnan(dndr)) = k; 
% Rx = -X;
% Ry = -Y;
% theta = acos((Rx.*V(1) + Ry.*V(2))./r);
% psi2 = ds.*sin(theta).*dndr./n;
% 
% figure(5)
% axis equal;
% % quiver(X,Y,Vkx,Vky)
% imagesc(psi2)
% title('deflection Meggit')
% 
% 
% 
% 
% 
% 
