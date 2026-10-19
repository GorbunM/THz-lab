A = crap;
nu = nuu(1:500);

%%

[peakVal, i0] = max(A);
nu0 = nu(i0);

halfVal = peakVal/2;
idx = find(A >= halfVal);

nuL = interp1(A(idx(1)-1:idx(1)), nu(idx(1)-1:idx(1)), halfVal);
nuR = interp1(A(idx(end):idx(end)+1), nu(idx(end):idx(end)+1), halfVal);

FWHM = nuR - nuL

%%

Q = nu0 / FWHM

%%

[pks, locs] = findpeaks(A, nu);

[mainPeak, imain] = max(pks);
mainNu = locs(imain);

sidePks = pks;
sidePks(imain) = [];

sideMax = max(sidePks);

SMSR_lin = mainPeak / sideMax;
SMSR_dB = 10*log10(SMSR_lin)

