close all

map = (tl.PP ~= 0);
p = [tl.PD(map) tl.PP(map)];
p = sortrows(p,1);

lensWidth = tl.centerLine - tl.leftEdge;
maxVal = max(p(:,2));
leftSpace = min(p(:,1));
rightSpace = lensWidth - max(p(:,1));

leftbox = [-lensWidth 0; -lensWidth maxVal;...
    (-lensWidth + rightSpace) maxVal; (-lensWidth + rightSpace) 0];
midbox = [-leftSpace 0; -leftSpace maxVal; leftSpace maxVal; leftSpace 0];
rightbox = [lensWidth 0; lensWidth maxVal;...
    (lensWidth-rightSpace) maxVal; (lensWidth - rightSpace) 0];
leftbox2 = [0 0; 0 maxVal; rightSpace maxVal; rightSpace 0];
rightbox2 = [lensWidth 0; lensWidth maxVal;...
    (lensWidth-leftSpace) maxVal; (lensWidth - leftSpace) 0];


p2 = p;
p2(:,1) = lensWidth-p2(:,1);
p3 = p2;
p3(:,1) = -p3(:,1);
pcomp = [p3; p2];

fp2 = fit(pcomp(:,1),pcomp(:,2),'poly2','Lower',[-inf -inf,1.64]);
fe2 = fit(p2(:,1),p2(:,2),'power2','Lower',[-inf -inf,1.64],'Upper',[0 inf inf]);
figure(2)
hold on
grid minor
x = linspace(0,lensWidth,1000);
yp2 = feval(fp2,x);
ye2 = feval(fe2,x(x>0));

plot(p2(:,1),p2(:,2),'*')
plot(x,yp2,'r','linewidth',1.4)
plot(x(x>0),ye2,'color','[0.3010 0.7450 0.9330]','linewidth',1.4)
plot([lensWidth;lensWidth],[0; maxVal],'k','linewidth',1)
plot([0;0],[0;maxVal],'k','linewidth',1)
xlabel('distance [pixels]')
ylabel('Phase-shift [lambda]')
legend('Data Points','f(x) = a*x^2 + b*x + c','f(x) = a*x^b+c')