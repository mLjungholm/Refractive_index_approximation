function AppPlotRefractiveGradient(imhandle)
% this.refractiveGradient = C.n;
% this.gradientD = C.r;
% C.plotNCurve
% C.plotN()
% if isempty(varargin)
%     scatter(this.r(2:end),this.n(2:end),'b*')
% else
%     scatter(this.r(varargin{1}+1),this.n(varargin{1}+1),'b*')
% end
% end
% function plotNCurve(this)
% for i = 1:this.nShells
%     x = [this.r(i); this.r(i+1)];
%     y = [this.n(i); this.n(i+1)];
%     plot(x,y,'r')
% end
n = this.refractiveGradient;
r = this.gradientD;
scatter(imhandle,r(2:end),n(2:end),'b*')
for i = 1:size(r,2)-1
    x = [this.r(i); this.r(i+1)];
    y = [this.n(i); this.n(i+1)];
    plot(imhandle,x,y,'r')
end

xmax = this.centerLine-this.leftEdge;
xt = linspace(0,xmax,10);
xt = round(xt);
imhandle.XTick = xt;

ymax = max(phase);
yt = linspace(0,ymax,10);
imhandle.YTick = yt;
end