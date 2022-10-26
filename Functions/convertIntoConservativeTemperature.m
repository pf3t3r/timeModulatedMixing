function [conservativeTemperature] = convertIntoConservativeTemperature(absoluteSalinity,temperature,p,time,label)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

filename = label;
conservativeTemperature = zeros(size(absoluteSalinity));

for i=1:length(time)
    conservativeTemperature(:,i) = gsw_CT_from_t(absoluteSalinity(:,i),temperature(:,i),p);
end
save(filename,'conservativeTemperature','temperature');

figure
plot(time,temperature(1,:),'DisplayName','T');
hold on
plot(time,conservativeTemperature(1,:),':','DisplayName','\Theta');
hold off
legend();
xlabel('time');
ylabel('T (C)');
title('Conservative temperature \Theta vs. temperature T');

end