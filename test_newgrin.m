close all

x = linspace(1,0,100);

n0 = 1.3;
n1 = 1.45;

[a,b,m] = getCoef(n0,n1);
modA = 2;
nFunc = @(x) (a*x.^2).*(1-sqrt(x)) + (b*x).*sqrt(x) + m ;
% nFunc2 = @(x) (a*x.^2).*(1-sqrt(x)).*(modA-(modA-1).*x) + (b*x).*sqrt(x) + m ;
nFunc2 = @(x) a.*x.^2 + m ;


figure(1)
hold on
% plot(x,nFunc(x),'b')
plot(x,nFunc2(x),'--r')

% figure(2)
% hold on
% plot(x,-x.^2,'b')
% % plot(x,sqrt(x).*(1.3-0.3.*x),'--r')

function [a,b,m] = getCoef(n0,n1)
a = (n0-n1);
b = n0-n1;
m = n1;
end
