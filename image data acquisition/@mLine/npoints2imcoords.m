function xy = npoints2imcoords(this)
if isempty(this.refractiveGradient)
%     x = nan; y = nan; n = nan;
xy = [];
    return
end
g_dist = this.centerLine - this.gradientD;
pnums = length(g_dist);
% closestP = zeros(pnums,1);
xy = zeros(pnums,2);
% y = zeros(pnums,1);
for pind = 1:pnums
    difd = abs(this.d - g_dist(pind));
    closestP = find(difd == min(difd));
    xy(pind) = this.coords(closestP,:);
%     x(pind) = this.coords(closestP,2);
end

end