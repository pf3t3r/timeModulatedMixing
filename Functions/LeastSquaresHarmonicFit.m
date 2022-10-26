function [A0, Amplitudes, Phases, X_mean, Fit, R2] = LeastSquaresHarmonicFit(t, f, X)

% CALCULATING THE LEAST SQUARES HARMONIC FIR TO A TIME SERIES.
%
% INPUT
%
% f     = Frequencies (rad/s) of the curves that re to be fitted.
% t     = Time (s);
% X     = Time series.

%
% OUTPUT
%
% Amplitudes;   Amplitudes related to frequencies.
% Phases;       Phases related to frequencies.
% X_mean;       The mean of the time series.
% LSHA;         The LEast squares Harmonic Fit.
%
% SCRIPT:
%
% We fit LSHA = A0 + Aa cos(wt) + Ab sin(wt).
% We first define the A and b matrices to construct A x = b and perform an 
% inversion to solve for Aa, Ab and A0 embedded in x.
% After that we use goniometric relationships to find:
% LSHA = Amplitudes * cos(wt - phases) for each frequency. 
% The sum over all frequencies is our fit.

%% Pre-calculations and defenitions:
% X_mean  = nanmean(X);       % Mean of time series.
X_mean = mean(X,'omitnan');
X0      = X - X_mean;       % Time series centered around 0.
N       = length(f);        % Number of frequencies to fit.
L       = length(t);
Ml      = 2*N + 1;          % Number of matrix elements.
A       = zeros(Ml,Ml);     % The A-matrix
b       = zeros(Ml,1);      % The b-matrix
f_mat   = f' * ones(1,L);   % Frequencies repeated for each time series.
t_mat   = ones(N,1) * t;    % Time series repeated for each Frequencie.
wt_mat  = f_mat .* t_mat;   % "omega t", i.e. the frequencie multiplied by time series.
Cn_mat  = cos(wt_mat);      % Cos(wt) for each frequency = [C1; C2; ...; Cn]
Sn_mat  = sin(wt_mat);      % Sin(wt) for each frequency = [S1; S2; ...; Sn]
CSn_mat = [Cn_mat;Sn_mat];  % = [C1; C2; ...; Cn; S1; S2; ...; Sn];

%% Filling the A matrix.
% To construct the middle part of A we use a matrix multiplication that 
% automatically includes a summation procedure due to the way matrices are 
% multiplied. This saves a for-loop:
A(1:2*N, 1:2*N)  = CSn_mat * CSn_mat' ;  % Part related to  related to (aj,bj) using [C1; C2; ...; Cn; S1; S2; ...; Sn] x [C1, C2, Cn, S1, S2, ..., Sn]. 
A(1:2*N, Ml)	 = sum(CSn_mat,2);       % Part related to A0.
A(Ml   , Ml)  	 = L;                    % Part related to A0. 
% WHY DOES THE ABOVE HAVE TO BE THE LENGTH OF THE TIME SERIES AND NOT JUST 1?????
% This was the only way it works, but now it works exactly.
A(Ml   ,1:2*N)   = sum(CSn_mat,2)';       % Create last row of A, i.e. with A0

%% Filling the b matrix.
X0_rep          = ones(2*N,1) * X0; % Repeat D0 time series.
% b(1:2*N,1)      = nansum( X0_rep .* CSn_mat, 2);
b(1:2*N,1)      = sum( X0_rep .* CSn_mat, 2,'omitnan');
clear X0_rep % The b-matrix related to (aj,bj).
% b(Ml,1)         = nansum(X0); % b-matric related to A0.
b(Ml,1) = sum(X0,'omitnan');

%% THE INVERSION
xn          = A \ b;
an          = xn(1:N);
bn          = xn(N+1 : 2*N );
A0          = xn(Ml);

%% CONSTRUCTING THE LSHA FIT:
Amplitudes      = sqrt(an.^2 + bn.^2);      % The amplitude for each cosine;         
Phases          = atan2(bn,an);             % HAS TO BE ATAN2, AND NOT JUST ATAN. The phases of the cosines.
Phases_mat      = Phases * ones(1,L);       % A matrix of phases.
Amplitudes_mat	= Amplitudes * ones(1,L);   % A marix of amplitudes.
Fit             = X_mean + A0 + sum( Amplitudes_mat .* cos(wt_mat - Phases_mat),1); % The LSHA fit.

%% ERROR ESTIMATE
SSE         = sum( (X - Fit).^2,'omitnan'); % Sum of squared errors
SST         = sum( (X - X_mean).^2,'omitnan'); % Sum of squared total variance
R2          = 1 - SSE/SST;             % Goodness of fit