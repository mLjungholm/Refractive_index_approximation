% Generates a interference pattern based on the calulated phase shift
% values of a ray source. 

% Input:
%       Source - traced source
%       Lambda - wavelength of source
%       r - radius of object
%       side - val = (-1,0,1). Determines what side of the centerpoint in
%       relation to the projection base vector. -1 for left side, 1 for
%       right side, 0 for both.

% Output:
%       c_phase_shift - caluclated phase shift in [m^-1]
%       c_phase_pos - caluclated position for phase hsift values along the
%       projection line
%       relative_phase_shift - value from [0,1] dependent on the
%       wavelength.
function [c_phase_shift, c_phase_pos, relative_phase_shift] = create_interference_pattern(source,lambda,r,side)
% remove all nan rays
key = ~isnan(source.projection(:,1));
c_phase_shift = source.phase(key);
c_phase_pos = source.projection(key,:);
% Calculate distance to center of projection line and what side of the
% center point each measured point lies.
for rayInd = 1:length(c_phase_pos(:,1))
    v = c_phase_pos(rayInd,:) - source.projection_center;
    c_phase_pos(rayInd,1) = norm(v);
    c_phase_pos(rayInd,2) = dot(v, source.projection_base_v);
end

% Sort out all values hot on the disiered side of the projection
switch side
    case 0
        c_phase_pos = c_phase_pos(:,1);
    case 1
        key = c_phase_pos(:,2) >= 0;
        c_phase_pos = c_phase_pos(key,1);
        c_phase_shift = c_phase_shift(key);
    case -1
        key = c_phase_pos(:,2) <= 0;
        c_phase_pos = c_phase_pos(key,1);
        c_phase_shift = c_phase_shift(key);
    otherwise
        fprintf('Error: not a valid side value. side must be -1,0 or 1. \n')
        return
end

% Sort positions and scale to radius.
comp = [c_phase_shift, c_phase_pos];
comp = sortrows(comp,2);
comp = comp.*r;
c_phase_shift = comp(:,1);
c_phase_pos = comp(:,2);

% Calculate relative phase shift to the refernece ray.
phaseFunc = @(x) (cos(x*2*pi/lambda + pi) + 1)/2;
relative_phase_shift = phaseFunc(c_phase_shift);

end