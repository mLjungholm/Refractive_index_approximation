function rayTrace2dGRIN_parallel(source,grin,stepSize,edge,n0)


arrayfun(@trace_ray,1:source.nRays);

    function trace_ray(rayInd)
        V = source.V(rayInd,:)';
        P = source.P(rayInd,:)';
                
        % Check if the ray hits the grin volume. Needs to be modified to
        % other shapes
        [p0,~,intersect] = circleIntersect(edge,P',V');
        if ~intersect
            source.terminated(rayInd) = 1;
            return
        end
        dt = stepSize;
        maxSteps = round(edge*2/dt*1.5);
        rayPath = zeros(maxSteps,2).*nan;
        diviation = zeros(maxSteps*2,1).*nan;
        stepError = zeros(maxSteps*2,1).*nan;
        rn = zeros(maxSteps*2,1).*nan;
        phase = 0;
        totalPath = 0;
        rayPath(1,:) = p0;
        stepNr = 1;
        P = p0';
        r0 = norm(p0);
        closest = findClosest(P);
        T = V.*grin.P(closest(2),closest(1));
        
        while ~source.terminated(rayInd) && stepNr < maxSteps
            closest = findClosest(P);
            rn(stepNr) = grin.P(closest(2),closest(1));
            dP = [grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            A = dt.*dP;
%             T = V.*grin.P(closest(2),closest(1));
            Bv = P + (dt/2).*T + 1/8*dt.*A;
            closest = findClosest(Bv);
            B = dt.*[grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            Cv = P + dt.*T + 1/2*dt.*B;
            closest = findClosest(Cv);
            C = dt.*[grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            R2 = P + dt.*(T + 1/6.*(A + 2.*B));
            T = T + 1/6.*(A + 4.*B + C);
            closest = findClosest(R2);
            V2 = T./grin.P(closest(2),closest(1));
            
            r1 = norm(R2);
            if r1 > edge
                source.terminated(rayInd) = 1;
            end
            ds = norm(P-R2);
            totalPath = totalPath + ds;
%             diviation(stepNr) = acos(dot(V',V2')./(norm(V)*norm(V2)));
            stepError(stepNr) = (ds-dt)/dt;
            phase = phase + ds*(rn(stepNr)-n0);
            rayPath(stepNr+1,:) = R2';
            P = R2;
%             V = V2;
%             r0 = r1;
            source.P(rayInd,:) = P;
            source.V(rayInd,:) = V;
            stepNr = stepNr + 1;
        end
        source.path{rayInd} = rayPath;
        source.phase(rayInd) = phase;
        source.diviation{rayInd} = diviation;
        source.totalPath(rayInd) = totalPath;
        source.stepError{rayInd} = stepError;
%         source.nPath{rayInd} = rn;
        source.P(rayInd,:) = P';
        source.V(rayInd,:) = V2';
        
        function closest = findClosest(P)
            [~,closestX] = min(abs(grin.X(1,:)-P(1)));
            [~,closestY] = min(abs(grin.Y(:,1)-P(2)));
            closest = [closestX,closestY];
        end
    end    
end