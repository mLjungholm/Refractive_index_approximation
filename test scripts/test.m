close all

na = 1.3;
nb = 1.4;
thetaA = 0:0.1:90;

figure(4)
thetaB = asind(na/nb*sind(thetaA));
psi = thetaA - thetaB;
plot(thetaA,psi)
