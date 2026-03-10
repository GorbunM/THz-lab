h  = 6.62607015e-34;      % J*s
kB = 1.380649e-23;        % J/K
c  = 299792458;           % m/s

T0 = 273;                 % K
T1 = 150;
T1 = T0 + T1;

T2 = 20;
T2 = T0 + T2; 

nu = linspace(0.1e12, 3e13, 1000);

Bnu1 = (2*h*nu.^3/c^2) ./ (exp(h*nu./(kB*T1)) - 1);
Bnu2 = (2*h*nu.^3/c^2) ./ (exp(h*nu./(kB*T2)) - 1);
dB = (Bnu1-Bnu2)*1e36;
itogo = dB .* A;


figure
hold on
% plot(nu*1e-12, Bnu1)
% plot(nu*1e-12, Bnu2)
plot(nu, itogo)