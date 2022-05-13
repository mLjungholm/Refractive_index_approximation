% Function normalizes all points to their own amplitude span.

function normalize_points(this)

for pInd = 1:this.pointNums
    pSpan = this.points(pInd,:);
    pSpan = pSpan - min(pSpan);
    this.points(pInd,:) = pSpan./max(pSpan);
end
this.normalized_points = 1;

this.smoothData
this.reset_line('peaks')

end