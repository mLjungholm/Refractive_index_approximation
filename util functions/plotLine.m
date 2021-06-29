function plotLine(p0,p1,color,figureNr)
figure(figureNr)
hold on 
plot([p0(1);p1(1)],[p0(2);p1(2)],color)
end