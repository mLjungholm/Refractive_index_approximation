function plotEllipse(ra,rb,rads,varargin)
if rads > 2*pi
    rads = 2*pi;
end
hold on
theta = linspace(-rads/2,rads/2,1000);
x = ra.*cos(theta);
y = rb.*sin(theta);
if isempty(varargin)
    plot(x,y,'linewidth',0.5,'Color',[0.6 0.6 0.6])
else
    plot(x,y,'linewidth',0.5,varargin{:})
end
end