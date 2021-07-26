function peaksDiff = get_phase_shift_from_peaks(peakVal,lambda)
halfsteps = (0:1:(length(peakVal)-2))./2.*lambda;
if peakVal(end-1) > 0.9
    p = halfsteps(end) + (1-peakVal(end))*lambda/2;
elseif peakVal(end-1) < 0.1
    p = halfsteps(end) + peakVal(end)*lambda/2;
else
    % If the value is here there is some sort of error
    disp('Error in getting the phase shift profile')
end
peaksDiff = [halfsteps';p];
if peakVal(1) > 0.1
    p0 = peakVal(1)*lambda/2;
    peaksDiff(1) = p0;
end
end