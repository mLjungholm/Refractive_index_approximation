function AppFitPhaseShiftProfile(this,fittype,restr)
p = this.PP(this.Pinc);
d = this.PD(this.Pinc);
w = this.centerLine - this.leftEdge;

% Mirror the data at the lens center
P = [p; flipud(p)];
dist_to_center = flipud(w-d);
D = [d; (flipud(d) +2.*dist_to_center)]; % Mirror around lens center
D = D - w; % shift the hwole set to be centered around 0

% f = fit(D,P,'poly2');
% x = linspace(-w,w,1000)';
% y = feval(f,x);
% miny = min(y(y>0));
% xpos = find(y(x>0) == miny) + find(x>0,1,'first') - 1;
% if  ~(miny == y(xpos))
%     this.fittype = 'failed fit';
%     return
% end
% this.edgeFit = x(xpos);
% this.fittype = 'poly2';

% switch fittype
%     case 'poly2'
if isequal(fittype,'poly2') || isequal(fittype,'poly4') || isequal(fittype,'poly6')
    if isequal(fittype,'poly2')
        if isequal(restr ,[-inf,inf])
            f = fit(D,P,'poly2');
        else
            f = fit(D,P,'poly2','Lower',[-inf -inf,restr(1)],'Upper',[inf inf,restr(2)]);
        end
    elseif isequal(fittype,'poly4')
        if isequal(restr ,[-inf,inf])
            f = fit(D,P,'poly4');
        else
            f = fit(D,P,'poly4','Lower',[-inf -inf,-inf,-inf,restr(1)],...
                'Upper',[inf inf,inf,inf,restr(2)]);
        end
    elseif isequal(fittype,'poly6')
        if isequal(restr ,[-inf,inf])
            f = fit(D,P,'poly6');
        else
            f = fit(D,P,'poly6','Lower',[-inf -inf,-inf,-inf,-inf,-inf,restr(1)],...
                'Upper',[inf inf,inf,inf,inf,inf,restr(2)]);
        end
    end
    x = linspace(-w,w,1000)';
    y = feval(f,x);
    miny = min(y(y>0));
    xpos = find(y(x>0) == miny) + find(x>0,1,'first') - 1;
    if  ~(miny == y(xpos))
        this.fittype = 'failed fit';
        return
    end
    this.fittype = fittype;
    this.edgeFit = x(xpos);
    sind = length(p) + 1;
    this.phasePoints = [P(sind:end),D(sind:end)];
elseif isequal(fittype,'power2')
    if restr == -inf
        f = fit(d,p,'power2');
    else
        f = fit(d,p,'power2','Lower',[-inf -inf,restr],'Upper',[0 inf inf]);
    end
    this.fittype = 'power2';
end
    this.phase_func = f;
end