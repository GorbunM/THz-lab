%% Interferogram
Npt = 301;
sigma  = linspace(400, 2000, Npt);  % см^-1

h   = 6.62607015e-27;      % erg·s
c   = 2.99792458e10;       % cm/s
kB  = 1.380649e-16;        % erg/K
c2  = h*c/kB;              % cm·K 
T   = 400;                 % K

Tsig = (2 * h * c^2) .* sigma.^3 ./ (exp(c2 * sigma / T) - 1);

halfway = 0.005;
Delta = linspace(-halfway, halfway, Npt);     % OPD, см

I = zeros(Npt,1);
for k = 1:numel(Delta)
    I(k) = trapz(sigma, Tsig.*(1+cos(2*pi*sigma*Delta(k))));
end

%% Noise
npt = 2001;
noiseLVL = 1;
noisefreq = 1000;

A = repmat(I, 1, npt);
% noise = rand(size(A)) - 0.5 * ones(size(A));    %flat
noise = [];                                       %low-freq 
for i = 1:Npt
    idx = [1, sort(randperm(npt-2,poissrnd(noisefreq))+1), npt];
    nos = nan(1, npt);
    nos(idx) = rand(size(idx))-0.5;
    nos = fillmissing(nos, 'linear');
    noise = [noise; nos];
end

Anoisy = A + 2 * max(A(:)) * noiseLVL * noise;

%% Average

Iaver = mean(Anoisy, 2);

%% Chopper
closed = 5;
open = 5;
chopper = [zeros(1, closed), ones(1, open)];
chopper = repmat(chopper, 1, ceil(npt/numel(chopper)));
chopper = chopper(1:npt);

Achop = A .* chopper + 2 * max(A(:)) * noiseLVL * noise;

%% Lockin 
t = linspace(0, npt, npt);
sinn = sin(2 * pi * t / (closed + open));
coss = cos(2 * pi * t / (closed + open));

f1 = Achop .* sinn;
g1 = Achop .* coss;

f2 = mean(f1, 2);
g2 = mean(g1, 2);

Ilockin = sqrt(f2 .^2 + g2 .^2);

%% Plots
figure
hold on
plot(Delta, I,"DisplayName", "original", "LineWidth", 1); 
plot(Delta, Iaver, "DisplayName", "averaging", "LineWidth", 1); 
% % plot(Delta, Ilockin, "DisplayName", "lock-in", "LineWidth", 1); 
plot(Delta, Ilockin./max(Ilockin) .*max(I), "DisplayName", "lock-in norm", "LineWidth", 1); 
% plot(Delta, Anoisy(:, randi(npt)), "DisplayName", "one-shot", "LineWidth", 1)
% plot(Delta, Anoisy(sub2ind([Npt npt], (1:Npt)', randi(npt, [Npt, 1]))));

leg = legend;
leg.FontSize = 24;
grid on;

xlabel('\Delta, cm'); ylabel('I');