function imCoords = npoints2imcoords(this) % [x,y]pairs
if isempty(this.refractiveGradient)
%     x = nan; y = nan; n = nan;
imCoords = [];
    return
end
gDist = this.centerLine - this.gradientD./this.pixelSize;
pnums = length(gDist);
% closestP = zeros(pnums,1);
imCoords = zeros(pnums,2);
% y = zeros(pnums,1);
lineDist = this.d;
if this.dataFliped
    lineDist = flipud(lineDist);
end
for pind = 1:pnums
    difd = abs(lineDist - gDist(pind));
    closestP = find(difd == min(difd));
    imCoords(pind,:) = this.coords(closestP,:);
%     x(pind) = this.coords(closestP,2);
end

end