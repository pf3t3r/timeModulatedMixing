function [SA_is,CT_is] = massStabilisation(SA_i,CT_i,p_i,time)
% massStabilisation
% This function takes as inputs interpolated conservative temperature CT_i,
% interpolated absolute salinity SA_i, and interpolated pressures p_i for
% one particular mooring and outputs the stabilised, interpolated
% conservative temperature CT_is and the stabilised, interpolated absolute
% salinity SA_is.
% Stabilisation is necessary to give non-negative buoyancies.

for i=1:length(time)
    [SA_is(:,i),CT_is(:,i)] = gsw_stabilise_SA_CT(SA_i(:,i),CT_i(:,i),p_i);
    disp(i);
end

end