function peaks2Phase(this)
if isempty(this.pks)
    fprintf('Error: No peaks have been selected for extraction. Use mLine.selectPeaks()')
    return
end

% phase = cell(this.imNums);
phase = [];
phaseD = pi;
layerPs = 2*pi/50;




pointNums = length(this.pks{1});
pd = this.pks{1};
ps = (0:1:pointNums-1)' .* phaseD;
f = fit(pd,ps,'poly2');
phase = [ps,pd];
for layerInd = 2:this.imNums
% for layerInd = 2:5
    pointNums = length(this.pks{layerInd});
    pd = this.pks{layerInd};
%     ps = ((0:1:pointNums-1)' .* phaseD) - layerPs*(layerInd-1);
ps = ((0:1:pointNums-1)' .* phaseD);
    pdelta = f(pd(1));
    ps = ps + pdelta;
    phase = [phase; ps, pd];
    phase = sortrows(phase,2);
    f = fit(phase(:,2),phase(:,1),'poly2');
end



h = figure(1);
hold on
grid on
% plot(pd,ps,'b.')
x = (0:this.centerLine/100:this.centerLine)';
plot(x,f(x),'r')
plot(phase(:,2),phase(:,1),'b.')
% ellipse_t = fit_ellipse( phase(:,2),phase(:,1),h);

end