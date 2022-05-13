function removeOutlier(this)

phaseInds = find(this.PP ~= 0);
p = [this.PD(phaseInds) this.PP(phaseInds)];

figure(1)
plot(p(:,1),p(:,2),'.')
grid minor
r = drawrectangle;
fprintf('Select exclution area. Press any key when done \n')
pause
r = r.Position;
close

span = [r(1) (r(1)+r(3))];
map1 = arrayfun(@(x) inSpan(x),p(:,1));
span = [r(2) r(2)+r(4)];
map2 = arrayfun(@(x) inSpan(x),p(:,2));
map = map1+map2;
map = map == 2;

this.PP(phaseInds(map)) = 0;

    function flag = inSpan(pi)
        if pi>= span(1) && pi <= span(2)
            flag = true;
        else
            flag = false;
        end
        fprintf('pd = %3.1f,  flag = %1u \n',pi,flag)
    end
end