function rayTrace2dGRIN_parallel(source,grin,stepSize,edge,n0)


arrayfun(@trace_ray,1:source.nRays);

    function trace_ray(rayInd)
        V = source.V(rayInd,:)';
        R = source.P(rayInd,:)';
                
        [p0,intersect] = circleIntersect(1,R',V');
        if ~intersect
            source.ended(rayInd) = 1;
            return
        end
        dt = stepSize;
        maxSteps = round(2/dt*1.5);
        rayPath = zeros(maxSteps,2).*nan;
        diviation = zeros(maxSteps*2,1).*nan;
        stepError = zeros(maxSteps*2,1).*nan;
        rn = zeros(maxSteps*2,1).*nan;
        phase = 0;
        totalPath = 0;
        rayPath(1,:) = p0;
        stepNr = 1;
        R = p0';
        r0 = norm(p0);
        closest = findClosest(R);
        T = V.*grin.P(closest(2),closest(1));
        
        while ~source.ended(rayInd) && stepNr < maxSteps
            closest = findClosest(R);
            rn(stepNr) = grin.P(closest(2),closest(1));
            dP = [grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            A = dt.*dP;
%             T = V.*grin.P(closest(2),closest(1));
            Bv = R + (dt/2).*T + 1/8*dt.*A;
            closest = findClosest(Bv);
            B = dt.*[grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            Cv = R + dt.*T + 1/2*dt.*B;
            closest = findClosest(Cv);
            C = dt.*[grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            R2 = R + dt.*(T + 1/6.*(A + 2.*B));
            T = T + 1/6.*(A + 4.*B + C);
            closest = findClosest(R2);
            V2 = T./grin.P(closest(2),closest(1));
            
            r1 = norm(R2);
            if r1 > edge
                source.ended(rayInd) = 1;
            end
            ds = norm(R-R2);
            totalPath = totalPath + ds;
%             diviation(stepNr) = acos(dot(V',V2')./(norm(V)*norm(V2)));
            stepError(stepNr) = (ds-dt)/dt;
            phase = phase + ds*(rn(stepNr)-n0);
            rayPath(stepNr+1,:) = R2';
            R = R2;
%             V = V2;
%             r0 = r1;
            source.P(rayInd,:) = R;
            source.V(rayInd,:) = V;
            stepNr = stepNr + 1;
        end
        source.PT{rayInd} = rayPath;
        source.phase(rayInd) = phase;
        source.diviation{rayInd} = diviation;
        source.totalPath(rayInd) = totalPath;
        source.stepError{rayInd} = stepError;
        source.nPath{rayInd} = rn;
        source.P(rayInd,:) = R';
        source.V(rayInd,:) = V2';
        
        function closest = findClosest(P)
            [~,closestX] = min(abs(grin.X(1,:)-P(1)));
            [~,closestY] = min(abs(grin.Y(:,1)-P(2)));
            closest = [closestX,closestY];
        end
    end    
end