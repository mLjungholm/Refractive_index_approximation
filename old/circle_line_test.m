close all
clear

c1 = 5;
l = [2 3 4 5 6];
ip = zeros(length(l),4);
d = zeros(length(l),1);

figure(1)
title('Ray Trace')
xlabel('[m]')
hold on; axis equal
viscircles([0,0],c1,'color','k','linewidth',0.5)
plot([-c1; c1].*1.2, [0;0],'color','k')
plot([0,0],[-c1; c1].*1.2,'color','k')

for i = 1:length(l)
    [dt, ipt] = intersect_dist(c1,l(i));
    if isnan(dt)
        continue;
    else
        ip(i,:) = ipt;
        d(i) = dt;
        plotLine(ipt);
    end
end

% circle: x^2 + y^2 = r^2
% verry simple equation since all the lines are vertical

function [d, ip] = intersect_dist(r,x)
if x >= r
    d = nan;
    ip = nan;
    return
end
% x is known, only y needs to be calculated
% y^2 = r^2 - x^2 => y = sqrt(r^2 - x^2);
% p1 = (x, y); p2 = (x, -y);
% d = ||p2-p1|| = ||(0,2*y)|| = sqrt((2*y)^2) = 2*y;
d = 2*sqrt(r^2 - x^2);

ip = [x sqrt(r^2 - x^2) x -sqrt(r^2 - x^2)];
end

function plotLine(ip)
plot(ip(1),ip(2),'ro')
plot(ip(3),ip(4),'ro')
plot([ip(1);ip(3)],[ip(2);ip(4)],'color','b')
end