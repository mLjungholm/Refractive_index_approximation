function plotCircle(r,rads,figureNr)
if rads > 2*pi
    rads = 2*pi;
end
figure(figureNr)
hold on
x = r.*cos(linspace(-rads/2,rads/2,1000));
y = r.*sin(linspace(-rads/2,rads/2,1000));
plot(x,y,'linewidth',0.5,'Color',[0.6 0.6 0.6])
end