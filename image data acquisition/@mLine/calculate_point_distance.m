%% Calculate distance from start of line
function calculate_point_distance(this)
this.d = zeros(this.pointNums,1);
p1 = this.coords(1,:);
% Deterimine the distance from start of the line to the mesurment point
for i = 1:this.pointNums
    p2 = this.coords(i,:);
    this.d(i) = norm(p2-p1);
end
end