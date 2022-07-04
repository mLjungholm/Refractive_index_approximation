function AppFlipData(this)
this.points = flipud(this.points);
this.gaussPoints = flipud(this.gaussPoints);
this.sgolayPoints = flipud(this.sgolayPoints);
this.d = this.d(end) - this.d;
this.d = flipud(this.d);
this.dataFliped = true;
this.reset_line('peaks')
this.centerLine = nan;
this.leftEdge = nan;

end