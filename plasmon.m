c = 29979245800; % speed of light in cm/s 

npt = 50;

e1 = 1;
e2 = 11.7;

omega = 2*pi*linspace(5e11, 6e12, npt);
KK = complex(zeros(4, numel(omega)));
RR = [280];

for idx = 1:numel(RR)
    Rsh = RR(idx)*(1e9/c^2);
    % Rsh = 300*(1e9/c^2);
    sigma0 = 1/Rsh;
    tau = 0.8e-13;
    % tau = RR(idx)*1e-15;
    sigma = sigma0./(1-1i.*omega*tau);
    
    K = complex(zeros(1, numel(omega)));
    Delta = zeros(1, numel(omega));
    K4 = complex(zeros(4, numel(omega)));
    Delta4 = zeros(4, numel(omega));
    idx4 = zeros(1, numel(omega));
    a4 = complex(zeros(4, numel(omega)));
    
    for i = 1:numel(omega)
    
        [K(:,i), Delta(:,i), K4(:,i), Delta4(:,i), idx4(:,i), a4(:,i)] = ...
            plasmonDispersionSolver(omega(i), sigma(i), e1, e2);
    
    end
    
    a = (K.*c./omega).^2;
    % a4 = (K4.*c./omega).^2;
    
    ro1 = e1./sqrt(e1-a);
    ro2 = e2./sqrt(e2-a);
    
    sig = (ro1 + ro2)*c/(4*pi);
    
    ro14 = e1./sqrt(e1-a4);
    ro24 = e2./sqrt(e2-a4);
    
    sig4 = (ro14 + ro24)*c/(4*pi);
    
    nu = omega*1e-12./(2*pi);
    KK(idx, :) = K;
    
end
%%
figure
hold on

for idx = 1:numel(RR)
    dispname = sprintf('\\tau = %g fs', RR(idx));
    % dispname = sprintf('R_{sh}=%g \\Omega /sq', RR(idx));
    colors = ["#08306B", "#08519C", "#0871B5", "#238B8D", "#41AB5D"];
    plot(real(KK(idx,:))*1e-4, nu, 'DisplayName', dispname)
    % plot(real(KK(idx,:))*1e-4./RR(idx), nu, 'DisplayName', dispname)
    % plot(real(KK(idx,:))*1e-4, nu, 'DisplayName', dispname, 'Color', colors(idx))
end

set(gca, 'LineWidth', 2, 'FontSize', 18)
set(findall(gcf, 'Type', 'Line'), 'LineWidth', 2)
ylabel('\nu [THz]')
% xlabel('k_p/R_{sh}')
xlabel('k_p')

xerr = 1e-2;
err = 1e-2;
% errorbar(exp(:,1), exp(:,2), xerr, xerr, err*1e-4, err*1e-4, 'ko', 'DisplayName', 'Experimental Data')
% plot(kIr/562.5, nuIr, '^')
% ylim([0, 2.5])
% xlim([0, 0.1])
% legend('Location','southeast')
grid on

%%

% writematrix(nuvsk', 'C:\Users\mgorbun\Documents\MATLAB\plasmonKvsNu.txt')