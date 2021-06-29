function rayTrace2dGRIN(source,grin,stepSize,halfStop)

for rayInd = 1:source.nRays
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
        continue
    end
    
    dt = stepSize;
    maxSteps = round(2/dt*1.5);
    stepNr = 3;
    rayPath = zeros(maxSteps,2);
    rayPath(:,1) = nan; rayPath(:,2) = nan;
    rayVpath = rayPath;
    raydV = rayPath;
    rayPath(1,:) = R';
    rayPath(2,:) = p;
    R = p';
    rayVpath(1,:) = V';
    N = -p./sqrt(p(1)^2+p(2)^2);
    [vNew, ~] = Snell(V', N, n1, n2);
    rayVpath(2,:) = vNew;
    V = vNew';
    inside = 0;
    while ~inside
        closest = findClosest(R);
        dP = [grin.dX(closest(2),closest(1));grin.dY(closest(2),closest(1))];
        if isnan(dP(1)) || isnan(dP(2))
            R = R + dt.*V;
            rayPath(stepNr,:) = R';
            rayVpath(stepNr,:) = V;
            stepNr = stepNr + 1;
        else
            inside = 1;
        end
    end
    
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
        if halfStop
            if R2(1) >= 0
                t = -R(1)/V2(1);
                rayPath(stepNr,:) = [0,R(2) + t*V2(2)];
                source.ended(rayInd) = 1;
                break
            end
        end
        if isnan(R2(1)) || isnan(R2(2))
            source.ended(rayInd) = 1;
            break
        end
        rayPath(stepNr,:) = R2';
        rayVpath(stepNr,:) = V';
        R = R2;
        V = V2;
        stepNr = stepNr + 1;
        if stepNr >= maxSteps
            source.ended(rayInd) = 1;
            break
        end
    end
    source.PT{rayInd} = rayPath;
end



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