function [P2, P1, f] = runDFT(noOfMeters,L,fs,rawU,rawV,FitU,FitV,resiU,resiV)
% runDFT: a function to quickly run a discrete Fourier transform (DFT)
% analysis of the data. The data is split up into raw data, tidal fit data,
% and residual data. We expect the DFT of the tidal fit to reproduce the
% tidal frequencies, while the residual data should reproduce some
% background signal that is due to mesoscale activity as well as background
% internal waves and incoherent features.
% P2: two-sided DFT (includes imaginary part)
% P1: one-sided DFT (imaginary part has been removed)
% f: frequency spectrum of the DFT
% noOfMeters: for how many vertical levels do we have data. This roughly
% corresponds to the number of instruments, but not quite since the
% instruments may move up or down during their time underwater.
% L: length (in time) of the dataset

rawU(isnan(rawU)) = 0;
rawV(isnan(rawV)) = 0;

dft_raw = zeros([2*noOfMeters,L]);
dft_fit = zeros([2*noOfMeters,L]);
dft_resi = zeros([2*noOfMeters,L]);

for i=1:noOfMeters
    dft_raw(i,:) = fft(rawU(i,:));
    dft_fit(i,:) = fft(FitU(i,:));
    dft_resi(i,:) = fft(resiU(i,:));
    dft_raw(i+noOfMeters,:) = fft(rawV(i,:));
    dft_fit(i+noOfMeters,:) = fft(FitV(i,:));
    dft_resi(i+noOfMeters,:) = fft(resiV(i,:));
end

P2_raw = abs(dft_raw/L);
P2_f = abs(dft_fit/L);
P2_res = abs(dft_resi/L);

P2 = [P2_raw; P2_f; P2_res];
P1 = P2(:,1:floor(L/2));       
f = (0:L/2-1)*fs/L;

end