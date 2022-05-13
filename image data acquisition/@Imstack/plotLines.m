function plotLines(this)
if isnan(this.lineNums)
    fprintf('ERROR: There are no data sampling lines. \n')
    return
end
figure('Name','Plot lines','NumberTitle','off')
C = this.imStack(:,:,1);
for lineInd = 1:size(this.mLines,1)
    if isempty(this.mLines{lineInd})
        continue
    end
    position = this.mLines{lineInd}.lineCoord(1,:);
    text_str = num2str(lineInd);
    C = insertText(C,position,text_str,'FontSize',20,'BoxColor',...
        'red','BoxOpacity',0.4,'TextColor','white');
end
imshow(C)
set(gcf,'position',this.figSize)
axis equal
if ~isempty(this.supLines) % Checks if there are any support lines and plots them
    for ind = 1:length(this.supLines)
        if ~isempty(this.supLines{ind})
            x = [this.supLines{ind}(1,1), this.supLines{ind}(2,1)];
            y = [this.supLines{ind}(1,2), this.supLines{ind}(2,2)];
        end
        line(x,y,'linewidth',1,'color','b')
    end
end
for lineInd = 1:size(this.mLines,1)
    if isempty(this.mLines{lineInd})
        continue
    end
    this.mLines{lineInd}.plotLine;
end
end