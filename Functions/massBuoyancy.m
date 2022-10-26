function [N2,pmid,z] = massBuoyancy(SA_is,CT_is,p_i,time,lat)
% massBuoyancy
% This function calculates the buoyancy frequency N^2 and pressure level
% that the buoyancy frequency is recorded at pmid. This is designed as an
% intermediate function to cut down on space in the main file.

N2 = [];
pmid = [];

for i=1:length(time)
    [x,y] = gsw_Nsquared(SA_is(:,i),CT_is(:,i),p_i');
    N2 = [N2 x];
    pmid = [pmid y];
    disp(i);
end

N2 = reshape(N2,[],length(time));
pmid = reshape(pmid,[],length(time));
z = gsw_z_from_p(pmid,lat);

end