% close all

figure(1)
hold on; grid on
s.plotTrace(1,'r')
s2.plotTrace(1,'b')

figure(2)
hold on
s1v = s.vPath{1}(:,1);
s2v = s2.vPath{1}(:,1);
s1v = s1v(~isnan(s1v));
s2v = s2v(~isnan(s2v));
plot(s1v,'r')
plot(s2v,'b')

figure(3)
hold on
s1n = s.nPath{1};
s2n = s2.nPath{1};
s1n = s1n(~isnan(s1n));
s2n = s2n(~isnan(s2n));
plot(s1n,'r')
plot(s2n,'b')
% figure(2)
% hold on
% grin.plot_nIndex_line(2)
% plot(linspace(1,0,100),nGradient(linspace(1,0,100)))