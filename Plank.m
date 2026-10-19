h  = 6.62607015e-34;      % J*s
kB = 1.380649e-23;        % J/K
c  = 299792458;           % m/s

T0 = 273;                 % K
T1 = 100;
T1 = T0 + T1;

T2 = 95;
T2 = T0 + T2; 

% nuu = linspace(1e12, 8e13, 1001);



B1 = (2*h*nu.^3/c^2) ./ (exp(h*nu./(kB*T1)) - 1);
B2 = (2*h*nu.^3/c^2) ./ (exp(h*nu./(kB*T2)) - 1);
% dB = (B1-B2)./(max(B1-B2));
dB = (B1-B2)*1e-2; % W/m^2 THz (0.01 sr)

% figure
% hold on
% plot(nu*1e-12, Bnu1)
% plot(nu*1e-12, Bnu2)
% plot(nu, dB)

%%
nuu = nu*1e-12;
lw = 2;

% nuuu = linspace(0.15, 20, 300);
% dB = interp1(nuu, dB, nuuu);
%%
ddB = dB.*Afig3;
% ddB = ddB(1:150);
figure
hold on
plot(nuu, ddB, 'DisplayName', '10 μm', 'LineWidth', lw)
plot(nuuu, dB, 'k:', 'DisplayName', 'Black Body', 'LineWidth', lw)
xlabel('\nu [THz]', 'FontSize', 14)
ylabel('L_\nu [W/(m^2\cdot THz)]', 'FontSize', 14)
ax = gca;
ax.FontSize = 24;
legend('FontSize', 24, 'Location', 'northwest')

box on
grid on
%%
nuu = nu(1:150);
ddB = ddB(1:150);
I = trapz(nuu, ddB) * pi * 0.05^2 * 1e3

