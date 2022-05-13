function plotPhase(this)

if nnz(this.PP) == 0
    fprintf('Error: No phase value has been extracted \n')
    return
end

% figure('Name','Phase plot','NumberTitle','off')
figure(1)
set(gcf,'position',[300,100,1000,800])
setFigure([0,this.d(end)])

for pInd = 1:this.pksNums
    if this.PP(pInd) ~= 0
%         txt = num2str(this.PL(pInd));
%         text(this.PD(pInd),this.PP(pInd),txt)
        plot(this.PD(pInd),this.PP(pInd),'b.','markerSize',15)
    end
end
maxX = this.centerLine;
maxY = max(this.PP);
startLine = [0 0; 0 maxY];
endLine = [maxX 0; maxX maxY];
plot(startLine(:,1),startLine(:,2),'k')
plot(endLine(:,1),endLine(:,2),'k')

if isnan(this.pixelSize)
    xtext = 'pixel distance';
else
    xtext = 'distance';
end
xlabel(xtext)
ylabel('Phase-shift [radians]')
end