function AppPlotPhaseGradient(this, imhandle)
colorlist = [0.3010 0.7450 0.9330;...
    0.4660 0.6740 0.1880;...
    0.4940 0.1840 0.5560;...
    0.8500 0.3250 0.0980;...
    0 0.4470 0.7410];
% d = this.PD(this.Pinc);
phase = this.PP(this.Pinc);

for groupInd = 1:length(this.S)
    key = this.Pinc & this.PS == groupInd;
    pd = this.PD(key);
    pp = this.PP(key);
    scatter(imhandle,pd,pp,'Marker','o','LineWidth',1,...
            'CData',colorlist(groupInd,:),'MarkerFaceColor','none')
%     scatter(imhandle,d,phase)
end

if ~isequal(this.fittype,'none')
    % Note that the fitted function is centered around 0. The shown values
    % in the plot needs to be mirrored and shifted
    x = linspace(0,this.edgeFit,1000);
    if isequal(this.fittype,'poly2')
        y = feval(this.phase_func,x);
    elseif isequal(this.fittype,'power2')
        y = feval(this.phase_func,x(x>0));
    end
    lw = this.centerLine-this.leftEdge;
    x = lw-x;
    plot(imhandle,x,y)
end

xmax = this.centerLine-this.leftEdge;
xt = linspace(0,xmax,10);
xt = round(xt);
imhandle.XTick = xt;

ymax = max(phase);
yt = linspace(0,ymax,10);
imhandle.YTick = yt;

% if ~isnan(this.centerZone) && ~isnan(this.centerLine)
%     v = [this.centerZone 0; this.centerZone ymax; this.centerLine ymax; this.centerLine 0];
%     f = [1 2 3 4];
%     patch('Faces',f,'Vertices',v,'FaceColor','[0.8500 0.3250 0.0980]','FaceAlpha',0.2,'Parent',imhandle);
% end

end