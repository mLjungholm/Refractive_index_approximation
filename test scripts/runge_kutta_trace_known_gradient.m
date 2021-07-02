function runge_kutta_trace_known_gradient(s,steps,gradient_type,n0,n1,r0)
ds = r0/steps;
switch gradient_type
    case 'parabolic'
        k = (n0-n1)/r0^2;
        n = @(r) k*(r(1)^2 + r(2)^2) + n1;
        D = @(r) 2*k*(k*(r(1)^2 + r(2)^2) + n1).*[r(1) r(2)];
%         D = @(r) 2*k.*[k.*(r(1)^3 + r(1)*r(2)^2) + n1*r(1);
%             k.*(r(2)^3 + r(2)*r(1)^2) + n1*r(2)];
    case 'linear'
        k = (n0-n1)/r0;
        n = @(r) k*sqrt(r(1)^2 + r(2)^2) + n1;
        D = @(r) k*(k + m/sqrt(r(1)^2 + r(2)^2)).*[r(1) r(2)];
        
end

arrayfun(@trace_ray,1:s.nRays);
    % Main tracing function
    function trace_ray(rayInd)
        % Initiate ray
        rayPath = zeros(steps*2,2).*nan;
        pathLength = 0;
        phaseShift = 0;
        p0 = s.P(rayInd,:);
        v0 = s.V(rayInd,:);
        % Check if volume intersection
        [p1,~,intersect] = circleIntersect(r0,p0,v0);        
        if ~intersect
            s.ended(rayInd) = 1;
            s.phase(rayInd) = 0;
            return
        end
        rayPath(1,:) = p1;
        p0 = p1;
        
        % Initiate tracing loop
        exitVolume = 0;
        rayStep = 2;
        while ~exitVolume && rayStep < steps*2
            [p1,v1] = rungekuttaTrace(p0,v0);
            r1 = sqrt(p1(1)^2 + p1(2)^2);
            if r1 > r0
                exitVolume = 1;
            end
            rayPath(rayStep,:) = p1;
            dt = norm(p0-p1);
            pathLength = pathLength + dt;
            phaseShift = phaseShift + dt*(n(p0)-n0);
            v0 = v1;
            p0 = p1;
            rayStep = rayStep + 1;
        end
        s.PT{rayInd} = rayPath;
        s.totalPath(rayInd) = pathLength;
        s.phase(rayInd) = phaseShift;
        s.P(rayInd,:) = p0;
        s.V(rayInd,:) = v0;
        s.ended(rayInd) = 1;
    end
        
% Runge-Kutta solution to the ray PDE (4th order) (O)~ds^5
    function [p1,v1] = rungekuttaTrace(p0,v0)
        A = ds.*D(p0);
        B = ds.*D(p0 + ds/2.*v0 + ds/8.*A);
        C = ds.*D(p0 + ds.*v0 + ds/2.*B);
        p1 = p0 + ds.*(v0 + 1/6.*(A + 2.*B));
        v1 = v0 + 1/6.*(A + 4.*B + C);
        v1 = v1./sqrt(v1(1)^2 + v1(2)^2);
    end
end