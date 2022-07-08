function AppCropStack(this,roi)
% roi = drawrectangle(imhandle);
% pause()
% [xmin, ymin, width, height]
this.cropRect = roi.Position;
% delete(roi)
imTemp = imcrop(this.imStackRaw(:,:,1), this.cropRect);
this.imSize = size(imTemp);
this.imStack = uint8(zeros(this.imSize(1),this.imSize(2),this.imNums));
for imInd = 1:this.imNums
    this.imStack(:,:,imInd) = imcrop(this.imStackRaw(:,:,imInd), this.cropRect);
end
end
