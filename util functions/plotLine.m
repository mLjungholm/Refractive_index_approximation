function plotLine(p0,p1,varargin)
if isempty(varargin)
    plot([p0(1);p1(1)],[p0(2);p1(2)],'Color',[0.6 0.6 0.6])
else
    plot([p0(1);p1(1)],[p0(2);p1(2)],varargin{:})
end
end