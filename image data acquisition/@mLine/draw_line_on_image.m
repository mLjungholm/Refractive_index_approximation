%% Plot the line on image
% Shift this to the imstack class
function draw_line_on_image(this)
x = [this.lineCoord(1,1), this.lineCoord(2,1)];
y = [this.lineCoord(1,2), this.lineCoord(2,2)];
line(x,y,'linewidth',1,'color','r')
end