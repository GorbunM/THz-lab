data = readmatrix("C:\\Users\\mgorbun\\Documents\\MATLAB\\FTIR\\STS_CNT2_2026.05.22\\STS_CNT2_2026.05.22\STS_CNT2_11Hz_12v_noflitr_raw.txt");

% data = readmatrix("C:\\Users\\mgorbun\\Documents\\MATLAB\\FTIR\\STS_CNT2_2026.05.22\\STS_CNT2_2026.05.22\yuki_sample_raw.txt");


x = data(:,1);   % mirror displacement
I = data(:,2);   % detector signal

delta = 2*x;     % OPD

I = I - mean(I); % remove DC

[~,i0] = max(abs(I));
I = circshift(I, round(length(I)/2)-i0);

% w = hann(length(I));
w = 1;
Iw = I .* w;

F = fftshift(fft(Iw));
S = abs(F);

d_delta = mean(diff(delta));
N = length(delta);

sigma = (-N/2:N/2-1)/(N*d_delta); % if delta in cm -> cm^-1

idx = sigma > 0;

plot(sigma(idx), S(idx))
xlabel('\sigma [cm^{-1}]')
ylabel('Amplitude')
grid on

















