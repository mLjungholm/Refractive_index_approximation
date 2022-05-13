fileID = fopen('incscape_selection2.svg');
text = textscan(fileID,'%s','Delimiter','');
text = text{1};
fileID = fclose(fileID);

match = regexp(text,'translate');
% [match, noMatch] = regexp(text,'translate','split','forceCellOutput');
cellID = find(~cellfun(@isempty,match));
coordShift = text{cellID}(match{cellID}+10:end-3);
coordShift = str2num(coordShift);

match = regexp(text,'^d="');
cellID = find(~cellfun(@isempty,match));
% Byt namn på edgeShift och edgePath. Värdena har hamnat fel!
edgePath = text{cellID(1)};
edgeShift = edgePath(regexp(edgePath,'c')+2:end-3);
edgeShift = str2num(edgeShift);
edgePath = edgePath(6:regexp(edgePath,'c')-2);
edgePath = str2num(edgePath);

centerLine = text{cellID(2)}(6:end-1);
centerLine = str2num(centerLine);

