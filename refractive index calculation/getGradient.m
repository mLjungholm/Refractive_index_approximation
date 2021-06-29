% Creates a gradient function from two given refractive indices and the
% raypath to create a gradient that gives the same raypath as a single
% refractive index.
close all

n0 = 1;
n1 = 1.2;
r0 = 3;
r1 = 2;
p0 = [r1,4];
v = [0,-1];

[ip0,ip1,intersect] = circleIntersect(r0,p0,v);
d = (ip1-ip0); d = sqrt(d(1)^2 + d(2)^2);

psi0 = d*(n1-n0);

figure(1)
hold on, axis equal, grid on
plotCircle(r0,1)
plotCircle(r1,1)
plotLine(p0,ip1,1)



function psi = linearStep(p0,v,d,dn,ds)
exitVol = 0;
maxStep = d/ds + 10;
rStep = 0;
psi = 0;
while ~exitVol && rStep <= maxStep
    
    
end
end