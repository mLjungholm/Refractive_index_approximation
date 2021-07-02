v0 = [0 -1];
p0 = [0 1.1];
rs = 1;
n0 = 1.3;
n1 = 1.4;
nProfile = 'linear';
nRays = 1;
steps = 10^5;
ds = rs/steps;
width = rs;

s_megg = source2d(p0,v0, nRays, width,n0,'half');
s_rk = source2d(p0,v0, nRays, width,n0,'half');
s_snell = source2d(p0,v0, nRays, width,n0,'half');


tic
meggit_trace(s_megg,ds,n0,n1,rs,nProfile,'meggit')
toc

tic
meggit_trace(s_rk,ds,n0,n1,rs,nProfile,'rk')
toc

tic 
meggit_trace(s_snell,ds,n0,n1,rs,nProfile,'snell')
toc

close all
figure(1)
hold on; axis equal; grid on
plotCircle(rs,2*pi,1)
plotLine([-rs 0], [rs, 0],'k',1)
plotLine([0,-rs], [0,rs],'k',1)
s_megg.plotTrace(1,'r')
s_rk.plotTrace(1,'b')
s_snell.plotTrace(1,'g')