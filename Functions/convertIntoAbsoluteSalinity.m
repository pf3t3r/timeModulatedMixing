function [absoluteSalinity] = convertIntoAbsoluteSalinity(practicalSalinity,time,p,lon,lat,label)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

filename = label;

absoluteSalinity = zeros(size(practicalSalinity));
% for j=1:length(noOfSensors)
    for i=1:length(time)
        absoluteSalinity(:,i) = gsw_SA_from_SP(practicalSalinity(:,i),p,lon,lat);
    end
% end
save(filename,'absoluteSalinity','practicalSalinity');

figure
plot(time,practicalSalinity(1,:),'DisplayName','S_P');
hold on
plot(time,absoluteSalinity(1,:),'DisplayName','S_A');
hold off
legend();
xlabel('Time');
ylabel('Salinity (g/kg)');
title('Absolute salinity S_A vs. practical salinity S_P');

end