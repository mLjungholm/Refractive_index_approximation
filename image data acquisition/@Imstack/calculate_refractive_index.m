function calculate_refractive_index(this,lineIndex)
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
% rightbox2 = rightbox + [1 0; 1 0; 1 0; 1 0].*rightSpace;
leftbox2 = [0 0; 0 maxVal; rightSpace maxVal; rightSpace 0];
rightbox2 = [lensWidth 0; lensWidth maxVal;...
    (lensWidth-leftSpace) maxVal; (lensWidth - leftSpace) 0];

f = [1 2 3 4];

outliers = [38; 64];
map = ones(size(p,1),1);
map(outliers) = 0;
p = p(map~= 0,:);

p2 = p;
p2(:,1) = lensWidth-p2(:,1);
p3 = p2;
p3(:,1) = -p3(:,1);
pcomp = [p3; p2];

f = fit(pcomp(:,1),pcomp(:,2),'poly2');
figure(2)
hold on
grid minor
x = linspace(0,lensWidth,1000);
y = feval(f,x);
plot(p2(:,1),p2(:,2),'.')
plot(x,y,'r')
plot([lensWidth;lensWidth],[0; maxVal],'k','linewidth',1)
plot([0;0],[0;maxVal],'k','linewidth',1)
xlabel('distance [pixels]')
ylabel('Phase-shift [lambda]')
end