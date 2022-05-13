function fit_phase_shift_profile(this)


map = (this.PP ~= 0);
p = [this.PD(map) this.PP(map)];
p = sortrows(p,1);

lensWidth = this.centerLine - this.leftEdge;
maxVal = max(p(:,2));

p2 = p;
p2(:,1) = lensWidth-p2(:,1);
p3 = p2;
p3(:,1) = -p3(:,1);
pcomp = [p3; p2];

fp2 = fit(pcomp(:,1),pcomp(:,2),'poly2');
fe2 = fit(p2(:,1),p2(:,2),'power2');
x = linspace(0,lensWidth,1000);
yp2 = feval(fp2,x);
ye2 = feval(fe2,x(x>0));


% lowerFit = -inf;
run = true;
while run
    plotPoly()
    fprintf('Do you whant to add f(0) lower restraint? \n')
    flag = yesno('');
    if ~flag
        fsave = input("Use polyFit or ExpoFit (f1 or f2)  Select[1/2]:");
        if fsave == 1
            this.phase_func = fp2;
        elseif fsave == 2
            this.phase_func = fe2;
        end
        run = false;
    else
        uval = input("Input lower bound: ");
        fitPoly(uval)
    end
    close(88)
end


    function fitPoly(uval)
        fp2 = fit(pcomp(:,1),pcomp(:,2),'poly2','Lower',[-inf -inf,uval]);
        fe2 = fit(p2(:,1),p2(:,2),'power2','Lower',[-inf -inf,uval],'Upper',[0 inf inf]);
        x = linspace(0,lensWidth,1000);
        yp2 = feval(fp2,x);
        ye2 = feval(fe2,x(x>0));
    end

    function plotPoly()
        figure(88)
        hold on
        grid minor
        plot(p2(:,1),p2(:,2),'*')
        plot(x,yp2,'r','linewidth',1.4)
        plot(x(x>0),ye2,'color','[0.3010 0.7450 0.9330]','linewidth',1.4)
        plot([lensWidth;lensWidth],[0; maxVal],'k','linewidth',1)
        plot([0;0],[0;maxVal],'k','linewidth',1)
        xlabel('distance [pixels]')
        ylabel('Phase-shift [lambda]')
        legend('Data Points','f(x) = a*x^2 + b*x + c','f(x) = a*x^b+c')
    end
end