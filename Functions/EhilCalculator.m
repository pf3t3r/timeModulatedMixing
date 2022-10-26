function [E_hil,E_hil_6dm,E_hil_Nbot2,E_hil_Nbot2_6dm,E_hil_U,E_hil_U_6dm,stdE,kurtE] = EhilCalculator(U,N2,N2_bot,power_hil)
% EhilCalculator
% Calculate the fucking Ehil for different moorings depending on the inputs
% provided.

% % Find density
% rho = gsw_rho(SA,CT,p);
% timeNum = meshgrid(datenum(time'),zq);
% clear CT SA;

% Apply floor to stratification
N2(N2 < 1e-11) = 1e-11;

% Set up the bottom buoyancy frequency
% This bottom buoyancy frequency estimate is based on the vertical average
% of the ten bottom-most estimates of the bouyancy frequency. Since these
% estimates are all themselves interpolated points this may not be the best
% method. An alternative method is to use the midpoint value of the two
% lowest sensors.

N2_bot(N2_bot<1e-9) = 1e-9;
N_bot = sqrt(N2_bot)';

% Normalise U
% U must be normalised such that the mean of U gives the single-valued
% parametrisation in space for E_hil
% Kunze (2017) states in Equation 6 that E_hil is a function of N_bot^2, U,
% and a roughness parameter h^2. Assuming that h^2 is contained within the
% existing parametrisation, we set the modulation according to the
% following:
alphaFactor = (N_bot.^2).*U;

% Normalise such that the mean of the alphaFactor = E_hil(x,y)
tuningFactor = 1./mean((alphaFactor - min(alphaFactor)) ./ (max(alphaFactor) - min(alphaFactor)));
normalisation = tuningFactor*(alphaFactor - min(alphaFactor)) ./ (max(alphaFactor) - min(alphaFactor));

% Normalise the individual contributions of Nbot^2 and U
normN = tuningFactor*(N_bot.^2 - min(N_bot.^2)) ./ (max(N_bot.^2) - min(N_bot.^2));
normU = tuningFactor*(U - min(U)) ./ (max(U) - min(U));

E_hil = normalisation*power_hil;
E_hil_Nbot2 = normN*power_hil;
E_hil_U = normU*power_hil;

E_hil_6dm = movmean(E_hil,149);
E_hil_Nbot2_6dm = movmean(E_hil_Nbot2,149);
E_hil_U_6dm = movmean(E_hil_U,149);

clear normN normU;

stdE = std(E_hil_6dm);
kurtE = kurtosis(E_hil_6dm);

end