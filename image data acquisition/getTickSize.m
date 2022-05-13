% Creates a ticksize depending on the axis scale and desired ticks

function  ticks = getTickSize(startVal, endVal, nrTicks)
tick = round((endVal-startVal)/nrTicks);
ndigits = numel(num2str(tick));
mods = [1 2 5 10].*10^(ndigits-1);
rest = abs(mods-tick);
ticks = mods(rest == min(rest));
end