function cropStack(this,imhandle)
fprintf('Select croping rectangle\n')
figure('Name','Select croping rectangle','NumberTitle','off')
imshow(this.imStackRaw(:,:,1),'Parent',imhandle);
set(gcf,'position',this.figSizeRaw)
axis equal
this.cropRect = getrect;
imTemp = imcrop(this.imStackRaw(:,:,1), this.cropRect);
this.imSize = size(imTemp);
this.imStack = uint8(zeros(this.imSize(1),this.imSize(2),this.imNums));
for imInd = 1:this.imNums
    this.imStack(:,:,imInd) = imcrop(this.imStackRaw(:,:,imInd), this.cropRect);
end
close();
ds = this.viewSize(3)/this.imSize(2);
if this.imSize(2)*ds > 800
    ds = this.viewSize(4)/this.imSize(1);
end
this.figSize = [this.viewSize(1) this.viewSize(2)...
    ds*this.imSize(2) ds*this.imSize(1)];
end

