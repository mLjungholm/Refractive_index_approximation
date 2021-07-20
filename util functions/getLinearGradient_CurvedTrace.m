% Simple function for finding the refractive index gradient in one shell
% that results in the same phase shift as the average refractive index "ns"
function n_test = getLinearGradient_CurvedTrace(r0,r1,n0,n1,phase,rayPath,varargin)
rayPath = rayPath(~isnan(rayPath(:,1)),:);
if phase <= 0
    n_test = n0;
    return
end
tresh = 0.01;
exit_flag = 0;
maxiter = 100;
iter = 1;
n_test = n1;
while exit_flag < tresh && iter < maxiter
    dphase = 0;
    k = (n0-n_test)/(r0-r1);
    nFunc = @(r) (r-r0)*k + n0;
    for rayInd = 1:length(rayPath(:,1))-1
        r = norm(rayPath(rayInd,:));
%         r = sqrt(rayPath(rayInd,1)^2 + rayPath(rayInd,2)^2);
        nr = nFunc(r);
        if isempty(varargin)
            ds = norm(rayPath(rayInd+1,:)-rayPath(rayInd,:));
            dphase = dphase + ds*(nr-n0);
        else
            dphase = dphase + varargin{1}*(nr-n0);
        end
    end
    acc = abs((dphase-phase)/phase);
    if acc > tresh
        if dphase < phase
            n_test = n_test+0.001;
        else
            n_test = n_test-0.001;
        end
    else
        exit_flag = 1;
    end
    iter = iter + 1;
end
end







