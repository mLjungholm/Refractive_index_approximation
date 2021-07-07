v0 = [0 -1];
p0 = [0 1.1];
rs = 1;
n0 = 1.3;
n1 = 1.5;
nProfile = 'parabolic';
% nProfile = 'linear';
% nProfile = 'eliptical';
k = (n1-n0)/rs;
n = @(r) k.*sqrt(rs^2 - r.^2) + n0;

% nProfile = 'luneburg';
nRays = 2;
steps = 10^4;
ds = rs/steps;
width = rs;
% half = 'half';
half = 0;

% g = GRIN2d(rs/10^3,nProfile,'matrix',n0,[n0 n1]);

s_megg = source2d(p0,v0, nRays, width,n0,half);
% s_megg_diff = source2d(p0,v0, nRays, width,n0,half);
s_rk = source2d(p0,v0, nRays, width,n0,half);
% s_grin = source2d(p0,v0, nRays, width,n0,half);


ray_trace(s_megg,ds,n0,n1,rs,nProfile,'meggit')

% ray_trace(s_megg_diff,ds,n0,n1,rs,nProfile,'megg_diff')

ray_trace(s_rk,ds,n0,n1,rs,nProfile,'rk')

% rayTrace2dGRIN_parallel(s_grin,g,ds,rs,n0);


close all
figure(1)
hold on; axis equal; grid on
plotCircle(rs,2*pi,1)
plotLine([-rs 0], [rs, 0],'k',1)
plotLine([0,-rs], [0,rs],'k',1)
s_megg.plotTrace(1,'r')
% s_megg_diff.plotTrace(1,'y')
s_rk.plotTrace(1,'b')
% s_grin.plotTrace(1,'b')

s_megg.getEndVals_single(1,'print');
s_rk.getEndVals_single(1,'print');


% figure(2)
% grid on
% plot(linspace(0,1,100),(n(linspace(0,1,100))),'b')