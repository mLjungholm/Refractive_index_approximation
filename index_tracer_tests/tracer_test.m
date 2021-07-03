v0 = [0 -1];
p0 = [0 1.1];
rs = 1;
n0 = 1.3;
n1 = 1.4;
% nProfile = 'parabolic';
% nProfile = 'linear';
% nProfile = 'eliptical';
nProfile = 'luneburg';
nRays = 2;
steps = 10^6;
ds = rs/steps;
width = rs;
% half = 'half';
half = 0;

% g = GRIN2d(rs/10^3,nProfile,'matrix',n0,[n0 n1]);

s_megg = source2d(p0,v0, nRays, width,n0,half);
% s_megg_diff = source2d(p0,v0, nRays, width,n0,'half');
s_rk = source2d(p0,v0, nRays, width,n0,half);
% s_snell = source2d(p0,v0, nRays, width,n0,'half');
% s_grin = source2d(p0,v0, nRays, width,n0,half);

tic
ray_trace(s_megg,ds,n0,n1,rs,nProfile,'meggit')
toc

% tic
% ray_trace(s_megg_diff,ds,n0,n1,rs,nProfile,'megg_diff')
% toc

tic
ray_trace(s_rk,ds,n0,n1,rs,nProfile,'rk')
toc

% tic 
% ray_trace(s_snell,ds,n0,n1,rs,nProfile,'snell')
% toc

% tic
% rayTrace2dGRIN_parallel(s_grin,g,ds,rs,n0);
% toc

close all
figure(1)
hold on; axis equal; grid on
plotCircle(rs,2*pi,1)
plotLine([-rs 0], [rs, 0],'k',1)
plotLine([0,-rs], [0,rs],'k',1)
s_megg.plotTrace(1,'r')
% s_megg_diff.plotTrace(1,'m')
s_rk.plotTrace(1,'b')
% s_snell.plotTrace(1,'g')
% s_grin.plotTrace(1,'k')

% figure(2)
% hold on; grid on
% plot(s_megg.diviation{1},'r')
% plot(s_megg_diff.diviation{1},'m')
% plot(s_rk.diviation{1},'b')
% plot(s_grin.diviation{1},'k')

% figure(3)
% hold on; grid on
% plot(s_megg.stepError{1},'r')
% plot(s_megg_diff.stepError{1},'m')
% plot(s_rk.stepError{1},'b')


