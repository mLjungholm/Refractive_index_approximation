function rayTrace2dGRIN_parallel(source,grin,stepSize,stopLine)


arrayfun(@trace_ray,1:source.nRays);

    function trace_ray(rayInd)
        %     % GRIN
        %     [X,Y,P] = create_2d_grin(0.01, 'linear', 'matrix');
        %     [px,py] = gradient(P,0.01);
        %     dX = P.*px;
        %     dY = P.*py;
        %
        %     % Ray
        %     V = [1; 0];
        %     R = [-2; 0.5];
        %     n1 = 1;
        %     n2 = 1;
        
        V = source.V(rayInd,:)';
        R = source.P(rayInd,:)';
        n1 = source.nStart;
        n2 = grin.nEdge;
        
        
        [p,intersect] = circleLineIntersect(1,R',V');
        
        if ~intersect
            source.ended(rayInd) = 1;
            return
        end
        
        dt = stepSize;
        maxSteps = round(2/dt*1.5);
        stepNr = 3;
        rayPath = zeros(maxSteps,2);
        rayPath(:,1) = nan; rayPath(:,2) = nan;
        vPath = rayPath;
        nPath = rayPath;
        rayVpath = rayPath;
        raydV = rayPath;
        rayPath(1,:) = R';
        rayPath(2,:) = p;
        vPath(1,:) = V;
        nPath(1) = n1;
        if p(2) == 0
           rayPath(3,:) = [p(1),stopLine];
           source.PT{rayInd} = rayPath;
           source.P(rayInd,:) = [p(1),stopLine];
           return
        end
        dp = p - R';
        source.OP(rayInd) = sqrt(dp(1)^2 + dp(2)^2)*n1;
        R = p';
        rayVpath(1,:) = V';
        N = -p./sqrt(p(1)^2+p(2)^2);
        [vNew, ~] = Snell(V', N, n1, n2);
        rayVpath(2,:) = vNew;
        V = vNew';
        vPath(2,:) = V;
        closest = findClosest(R);
        nPath(2) = grin.P(closest(2),closest(1));
%         inside = 0;
%         extraStep = 1;
%         while ~inside
%             closest = findClosest(R);
%             dP = [grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
%             if isnan(dP(1)) || isnan(dP(2))
%                 R = R + dt.*V;
%                 rayPath(stepNr,:) = R';
%                 rayVpath(stepNr,:) = V;
%                 stepNr = stepNr + 1;
%             else
%                 inside = 1;
%             end
%             extraStep = extraStep + 1;
%             if extraStep > 4
%                 rayPath = [];
%                 source.ended(rayInd) = 1;
%                 return
%             end
%         end
        
        while ~source.ended(rayInd)
            closest = findClosest(R);
            dP = [grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            raydV(stepNr,:) = dP';
            A = dt*dP;
            Bv = R + (dt/2)*V + 1/8*dt*A;
            closest = findClosest(Bv);
            B = dt*[grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            Cv = R + dt*V + 1/2*dt*B;
            closest = findClosest(Cv);
            C = dt*[grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
            R2 = R + dt*(V + 1/6*(A + 2*B));
            V2 = V + 1/6*(A + 4*B + C);
            closest = findClosest(R);
            if R2(2) <= stopLine
                t = (stopLine-R(2))/V2(2);
                R2 = [R(1) + t*V2(1),stopLine]';
                dp = R2 - R;
                source.OP(rayInd) = source.OP(rayInd) + sqrt(dp(1)^2 + dp(2)^2)*grin.P(closest(2),closest(1));
                rayPath(stepNr,:) = R2';
                source.ended(rayInd) = 1;
                break
            end
            if isnan(R2(1)) || isnan(R2(2))
                t = (stopLine-R(2))/V(2);
                R2 = [R(1) + t*V(1),stopLine]';
                dp = R2 - R;
                source.OP(rayInd) = source.OP(rayInd) + sqrt(dp(1)^2 + dp(2)^2)*grin.P(closest(2),closest(1));
                rayPath(stepNr,:) = R2';
                source.ended(rayInd) = 1;
                break
            end
            dp = R2 - R;
            source.OP(rayInd) = source.OP(rayInd) + sqrt(dp(1)^2 + dp(2)^2)*grin.P(closest(2),closest(1));
            rayPath(stepNr,:) = R2';
            rayVpath(stepNr,:) = V';
            vPath(stepNr,:) = V';
            nPath(stepNr) = grin.P(closest(2),closest(1));
            R = R2;
            V = V2;
            source.P(rayInd,:) = R;
            source.V(rayInd,:) = V;
            stepNr = stepNr + 1;
            if stepNr >= maxSteps
                source.ended(rayInd) = 1;
                break
            end
        end
        source.PT{rayInd} = rayPath;
        source.vPath{rayInd} = vPath;
        source.nPath{rayInd} = nPath;
        function closest = findClosest(P)
            [~,closestX] = min(abs(grin.X(1,:)-P(1)));
            [~,closestY] = min(abs(grin.Y(:,1)-P(2)));
            closest = [closestX,closestY];
        end
        
        function [p,intersect] = circleLineIntersect(r,p0,v)
            L = -p0;
            tca = v(1)*L(1) + v(2)*L(2);
            if tca < 0
                intersect = 0;
                p = nan;
                return
            end
            d = sqrt(L(1)^2 + L(2)^2 - tca^2);
            if d < 0 || d > r
                intersect = 0;
                p = nan;
                return
            end
            thc = sqrt(r^2-d^2);
            t0 = tca - thc;
            p = p0 + t0*v;
            intersect = 1;
        end
    end



    
end