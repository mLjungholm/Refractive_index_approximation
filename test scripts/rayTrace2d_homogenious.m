function rayTrace2d_homogenious(source)


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
        n2 = 1.5;
        
        
        [p,intersect] = circleLineIntersect(1,R',V');
        
        if ~intersect
            source.ended(rayInd) = 1;
            return
        end
        
        maxSteps = 4;
        rayPath = zeros(maxSteps,2);
        rayPath(:,1) = nan; rayPath(:,2) = nan;
        rayVpath = rayPath;
        rayPath(1,:) = R';
        rayPath(2,:) = p;
        dp = p - R';
        source.OP(rayInd) = sqrt(dp(1)^2 + dp(2)^2)*n1;
        R = p';
        rayVpath(1,:) = V';
        N = -p./sqrt(p(1)^2+p(2)^2);
        [vNew, ~] = Snell(V', N, n1, n2);
        rayVpath(2,:) = vNew;
        V = vNew';
       
        t = -R(1)/V(1);
        y = R(2) + t*V(2);
        R2 = [0,y]';
        rayPath(3,:) = R2';
        dp = R - R2;
        source.PT{rayInd} = rayPath;
        source.OP(rayInd) = source.OP(rayInd) + sqrt(dp(1)^2 + dp(2)^2)*n2;
        source.P(rayInd,:) = R2';
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