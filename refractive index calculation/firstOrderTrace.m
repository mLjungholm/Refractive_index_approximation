% Initial approximation of the refractive index in each shell of the slice
% s - ray source
% pp - phase shift positions. order [edge-center]
% ps - phase shifts
% n0 - refractive index of surounding medium
function firstOrderTrace(s,g)
shellN = g.nShells;
n = g.n;
rayInd = 1;

n(rayInd+1) = firstTrace();
for i = 2:shellN-1
    rayInd = rayInd + 1;
    n(rayInd+1) = secondTrace();
end
%     Approximation for outmost ray. i.e the first shell
    function nR = firstTrace()
        % Outer shell intersection points
        [ip0,ip1,intersect] = circleIntersect(g.shellR(rayInd),s.P(rayInd,:),s.V(rayInd,:));
%         plotLineIntersect(s.P(rayInd,:),ip1)
        if ~intersect
            disp('Error: no sircle intersection')
            return
        end
        if ip1 == ip0
            disp('No intersect: ray tangent')
            nR = g.n0;
            return
        end
        dv = (ip1-ip0);
        d = sqrt(dv(1)^2+dv(2)^2);
        nR = g.n0 + (s.shellR(rayInd)/d); % Calculate refractive index from known phase shift and distance
    end

    % Approximation of consecutive shells
    function nR = secondTrace()
        nShells = rayInd; % Number of shells to be taken into acount
        dShells = zeros(nShells,1); % Distance traveled in each shell.
        for shellInd = 1:nShells
            try
                [ip0,ip1,~] = circleIntersect(pp(shellInd),s.P(rayInd,:),s.V(rayInd,:));
            catch
                disp(rayInd)
                testflag = 1;
            end
%             plotLineIntersect(s.P(rayInd,:),ip0)
%             plotLineIntersect(s.P(rayInd,:),ip1)
            dv = (ip1-ip0);
            d = sqrt(dv(1)^2+dv(2)^2); % distance traveled in shell
            dShells(shellInd) = d;
        end
        % Substract the distance from inner shells.
        for shellInd = 1:nShells-1
            dShells(shellInd) = dShells(shellInd) - dShells(shellInd+1); 
        end
        
        dPhase = zeros(nShells,1);
%         dPhase(1) = dShells(1) * (n(2)-n0);
        for shellInd = 1:nShells-1
            dPhase(shellInd) = dShells(shellInd)*(n(1+shellInd) - n(1));
        end
        dPhase(rayInd) = ps(rayInd+1) - sum(dPhase(1:rayInd-1));
        nR = dPhase(end)/dShells(end) + n0;
    end

    function plotLineIntersect(p,ip1)
%         plot(ip0(1),ip0(2),'o')
%         plot(ip1(1),ip1(2),'o')
        plot([p(1);ip1(1)],[p(2);ip1(2)],'r')
    end
end






