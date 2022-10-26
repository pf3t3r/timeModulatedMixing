function [A0,Amps,Phases,Xmean,Fit,R2,resi] = runLSHA(time,wn,u,v,noOfMeters)
% runLSHA: quickly run the LSHA function on u and v according to the number
% of different current-meters / ADCPs and automatically output a file
% containing the residual.

u(isnan(u)) = 0;
v(isnan(v)) = 0;

for i=1:noOfMeters
    [A0.u(i), Amps.u(i,:), Phases.u(i,:), Xmean.u(i), Fit.u(i,:), R2.u(i)] = LeastSquaresHarmonicFit(time, wn, u(i,:));
    [A0.v(i), Amps.v(i,:), Phases.v(i,:), Xmean.v(i), Fit.v(i,:), R2.v(i)] = LeastSquaresHarmonicFit(time, wn, v(i,:));
end

resi.v = v - Fit.v;
resi.u = u - Fit.u;
resi.v(isnan(resi.v)) = 0;  % must set NaN = 0 for the filtering to work
resi.u(isnan(resi.u)) = 0;

end