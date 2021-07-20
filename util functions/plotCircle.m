function plotCircle(r,rads,varargin)
if rads > 2*pi
    rads = 2*pi;
end
hold on
x = r.*cos(linspace(-rads/2,rads/2,1000));
y = r.*sin(linspace(-rads/2,rads/2,1000));
if isempty(varargin)
    plot(x,y,'linewidth',0.5,'Color',[0.6 0.6 0.6])
else
    plot(x,y,'linewidth',0.5,'Color',varargin{1})
end
end