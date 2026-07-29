%% 1. Constants (CGS units)

c = 29979245800; % speed of light in cm/s 

%% 2. Frequency range

npt = 300;
nu = linspace(1.5e11, 40e12, npt);  % frequency in Hz

omega = 2*pi*nu;  % angular frequency in rad/s

%% 3. Layer definitions via Lorentz/Drude parameters

% Layer types: -1 -- given
%               0 -- air/vacuum
%               1 -- Drude
%               2 -- Lorenz
%               3 -- Drude + Lorenz 

% Layer templates:
   
    % Given = struct('type', -1, ...
    %              'd', 2e-3, ...
    %              'n0', 0, ...
    %              'eps1', [eps1A, eps1B], ...
    %              'eps2', [eps2A, eps2B], ...
    %              'omega', omega);

    % Air = struct('type', 0, ...
    %              'd', 1, ...
    %              'n0', 1);

    % DrudeLayer = struct('type', 1, ...
    %                     'd', 1e-4, ...
    %                     'n0', 1, ...
    %                     'omega_p', 2*pi*2e14, ...
    %                     'tauD', 2*pi*1e13);
    
    % LorentzLayer = struct('type', 2, ...
    %                       'd', 5e-5, ...
    %                       'n0', 1, ...
    %                       'tau', [2*pi*1e13, 2*pi*5e12], ...
    %                       'omega_0', [2*pi*1.5e14, 2*pi*1.2e14], ...
    %                       'magn', [10, 1]);

    % DrudeLorentzLayer = struct('type', 3, ...
    %                         'd', 1e-4, ...
    %                         'n0', sqrt(3.85), ...
    %                         'omega_p', 4.5e13, ...
    %                         'tauD', 1 / 3.9e12, ...
    %                         'tau', 1 / 3.9e12, ...
    %                         'omega_0', 7.55e13, ...
    %                         'magn', 1);

p = 3.5/3.5;
m = 10;

    Air = struct('type', 0, ...
             'd', 1, ...
             'n0', 1);

    Di = struct('type', 1, ...
                'd', 10e-4*p, ...
                'n0', 3.5/p, ...
                'omega_p', 0, ...
                'tauD', 0);
    
    Au = struct('type', 1, ...
                'd', 2e-5, ...
                'n0', 1, ...
                'omega_p', 2*pi*2.15e15, ...
                'tauD', 1e-14);
    
    TL = struct('type', 2, ...
                'd', 1e-5, ...
                'n0', 1, ...
                'tau', 1e-12/(2*pi), ...
                'omega_0', 2*pi*10.7e12, ... %2*pi*1.05e13, ...
                'magn', 4);

%% 4. Stack construction

structure = {Air, TL, Di, Au, Air};
N = length(structure);

%% 5. Compute complex refractive index for each layer at all frequencies

n_layers = zeros(N, npt);
eps_layers = zeros(N, npt);
for j = 1:N
    curlayer = structure{j};
    eps = curlayer.n0 ^ 2;
    if curlayer.type == -1
        eps = curlayer.eps1 + 1j*curlayer.eps2;                             % Given
    
    elseif curlayer.type == 1
        w_p = curlayer.omega_p;
        g   = 1 / curlayer.tauD;
        eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);                    % Drude model
    
    elseif curlayer.type == 2
        for peak = 1:length(curlayer.omega_0)
            g   = 1 / curlayer.tau(peak);
            w_0 = curlayer.omega_0(peak);
            A_0 = curlayer.magn(peak);
            eps = eps + A_0 * (w_0^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Lorentz model
        end
    
    elseif curlayer.type == 3
        w_p = curlayer.omega_p;
        A_0 = curlayer.magn;
        g   = 1 / curlayer.tauD;
        eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);
    
        for peak = 1:length(curlayer.omega_0)
            g   = 1 / curlayer.tau(peak);
            w_0 = curlayer.omega_0(peak);
            eps = eps + A_0(peak) * (w_0^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Drude-Lorentz model
        end
        
        n_layers(j, :) = sqrt(eps);  % complex refractive index
        eps_layers(j, :) = eps;  % complex 
    
    end
    n_layers(j, :) = sqrt(eps);  % complex refractive index
    eps_layers(j, :) = eps;  % complex refractive index
end

%% 6. Assemble transfer matrices


T = repmat(eye(2), 1, 1, npt);  % start with 2×2×N identity matrices

for j = 2:N-1
    curlayer = structure{j};
    d = curlayer.d;
    nL = n_layers(j-1, :);
    nR = n_layers(j,   :);

    % Interface matrix M
    M = zeros(2,2,npt);
    M(1,1,:) = 1 + nL ./ nR;
    M(1,2,:) = 1 - nL ./ nR;
    M(2,1,:) = 1 - nL ./ nR;
    M(2,2,:) = 1 + nL ./ nR;
    M = 0.5 * M;
    
    % Multiply: T = T * B * P
    T = pagemtimes(T, M);

    delta = omega .* nR .* d / c;  % 1×N

    % Propagation matrix P
    P = zeros(2,2,npt);

    P(1,1,:) = exp(1i*delta);
    P(2,2,:) = exp(-1i*delta);
    
    T = pagemtimes(T, P);
end

% Final interface (N−1 to N)
nL = n_layers(N-1,:);
nR = n_layers(N,  :);

M = zeros(2,2,npt);
M(1,1,:) = 1 + nL ./ nR;
M(1,2,:) = 1 - nL ./ nR;
M(2,1,:) = M(1,2,:);
M(2,2,:) = M(1,1,:);
M = 0.5 * M;

T = pagemtimes(T, M);  % final multiply

%% 7. Compute reflection and transmission spectra 
Ein = squeeze(T(:,2,:));  % incoming field from right-side excitation

refl = abs(Ein(1,:) ./ Ein(2,:)).^2;
trans = abs(1 ./ Ein(2,:)).^2;
phase = angle(Ein(1,:) ./ Ein(2,:));

A = 1 - refl - trans;  % absorption (if any)

%% Plot

nuu = nu*1e-12;

figure
hold on
plot(nuu, A, 'k', 'DisplayName','Absorption', 'LineWidth', 2)
% plot(nuu, refl, 'r', 'DisplayName','Reflectance','LineWidth', 2)
% plot(nuu, trans, 'b', 'DisplayName','T','LineWidth', 2)
ylim([0, 1])

% plot(nuu, phase, 'DisplayName','$\phi$')
% ylim([0, pi/2])

% plot(nuu, real(eps_layers(2,:)), 'DisplayName', '$\varepsilon^{\prime}$')
% plot(nuu, imag(eps_layers(2,:)), 'DisplayName', '$\varepsilon^{\prime\prime}$')

xlabel('\nu, THz')
ylabel('Emissivity')
ax = gca;
ax.FontSize = 24;

% legend('Interpreter','latex', 'FontSize', 24, 'Location','east')

box on
grid on

%% 8. Field and loss profiles from TMM

% Frequencies for comparison
[~, k_res] = max(A);              % resonance frequency
[~, k_off] = min(abs(nu - 8e12)); % off-resonance reference, change if needed

k_list = [k_res, k_off];
labels = {'resonance', 'off-resonance'};

profiles = struct([]);

for kk = 1:length(k_list)

    k_ind = k_list(kk);

    % Complex reflection amplitude, not reflectance
    r = Ein(1,k_ind) / Ein(2,k_ind);

    % Initial amplitudes in incident medium:
    % E(z) = E_plus exp(i beta z) + E_minus exp(-i beta z)
    amp = [1; r];

    z_all = [];
    E2_all = [];
    ploss_all = [];
    layer_id = [];

    z_offset = 0;

    for j = 2:N-1

        % Interface from layer j-1 to layer j
        nL = n_layers(j-1,k_ind);
        nR = n_layers(j,k_ind);

        M = 0.5 * [1 + nL/nR, 1 - nL/nR; ...
                   1 - nL/nR, 1 + nL/nR];

        % Amplitudes just inside layer j
        amp = M \ amp;

        % Local coordinate inside layer j
        d = structure{j}.d;

        if j == 2 || j == 4
            Nz = 300;     % thin TL / Au
        else
            Nz = 1000;    % dielectric spacer
        end

        zloc = linspace(0, d, Nz);

        beta = omega(k_ind) * nR / c;

        E = amp(1) * exp(1i * beta * zloc) + ...
            amp(2) * exp(-1i * beta * zloc);

        E2 = abs(E).^2;

        % CGS loss density up to physical normalization:
        % q(z) = omega/(8*pi) * Im(eps) * |E|^2
        eps_im = imag(eps_layers(j,k_ind));
        ploss = omega(k_ind) / (8*pi) * eps_im * E2;

        z_all = [z_all, z_offset + zloc];
        E2_all = [E2_all, E2];
        ploss_all = [ploss_all, ploss];
        layer_id = [layer_id, j * ones(size(zloc))];

        % Propagate amplitudes to the right boundary of layer j
        P = [exp(1i * beta * d), 0; ...
             0, exp(-1i * beta * d)];

        amp = P * amp;

        z_offset = z_offset + d;
    end

    profiles(kk).k_ind = k_ind;
    profiles(kk).nu_THz = nu(k_ind) * 1e-12;
    profiles(kk).z = z_all;
    profiles(kk).E2 = E2_all;
    profiles(kk).ploss = ploss_all;
    profiles(kk).layer_id = layer_id;
end

%% 9. Plot |E|^2 profile: resonance vs off-resonance

figure
hold on

plot(profiles(1).z * 1e4, profiles(1).E2, 'k', ...
     'LineWidth', 2, ...
     'DisplayName', ['resonance, ', num2str(profiles(1).nu_THz, '%.2f'), ' THz'])

plot(profiles(2).z * 1e4, profiles(2).E2, '--k', ...
     'LineWidth', 2, ...
     'DisplayName', ['off-resonance, ', num2str(profiles(2).nu_THz, '%.2f'), ' THz'])

z_TL = structure{2}.d * 1e4;
z_Di = (structure{2}.d + structure{3}.d) * 1e4;

xline(z_TL, ':', 'TL/Di', 'LineWidth', 1.5, 'LabelVerticalAlignment','top')
xline(z_Di, ':', 'Di/Au', 'LineWidth', 1.5, 'LabelVerticalAlignment','top')

xlabel('z, \mum')
ylabel('|E|^2 / |E_0|^2')
title('Electric-field intensity profile')

% legend('Location','best')
ax = gca;
ax.FontSize = 24;
box on
grid on

%% 10. Zoom: resonant layer only

figure
hold on

idx_TL_res = profiles(1).layer_id == 2;
idx_TL_off = profiles(2).layer_id == 2;

plot(profiles(1).z(idx_TL_res) * 1e4, profiles(1).E2(idx_TL_res), ...
     'k', 'LineWidth', 2, ...
     'DisplayName', ['resonance, ', num2str(profiles(1).nu_THz, '%.2f'), ' THz'])

plot(profiles(2).z(idx_TL_off) * 1e4, profiles(2).E2(idx_TL_off), ...
     '--k', 'LineWidth', 2, ...
     'DisplayName', ['off-resonance, ', num2str(profiles(2).nu_THz, '%.2f'), ' THz'])

xlabel('z, \mum')
ylabel('|E|^2 / |E_0|^2')
title('Electric-field intensity inside resonant layer')

% legend('Location','best')
ax = gca;
ax.FontSize = 24;
box on
grid on

%% 11. Power-loss density profile at resonance

figure
plot(profiles(1).z * 1e4, profiles(1).ploss, 'k', 'LineWidth', 2)
hold on

xline(z_TL, ':', 'TL/Di', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom')
xline(z_Di, ':', 'Di/Au', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom')

xlabel('z, \mum')
ylabel('q(z), arb. units')
title(['Power-loss density at \nu = ', num2str(profiles(1).nu_THz, '%.2f'), ' THz'])

ax = gca;
ax.FontSize = 24;
box on
grid on

%% 12. Integrated absorption/loss by layers at resonance

layers_to_integrate = [2, 3, 4];
layer_names = {'TL', 'Dielectric', 'Au'};

loss_int = zeros(size(layers_to_integrate));

for jj = 1:length(layers_to_integrate)

    j = layers_to_integrate(jj);
    idx = profiles(1).layer_id == j;

    z_j = profiles(1).z(idx);
    q_j = profiles(1).ploss(idx);

    loss_int(jj) = trapz(z_j, q_j);
end

loss_frac = loss_int / sum(loss_int);

figure
bar(loss_frac)

set(gca, 'XTick', 1:length(layer_names), 'XTickLabel', layer_names)

ylabel('Fraction of absorbed power')
title(['Absorption distribution at \nu = ', num2str(profiles(1).nu_THz, '%.2f'), ' THz'])

ax = gca;
ax.FontSize = 24;
box on
grid on

%% 13. Print numbers

fprintf('\nField/loss analysis at resonance:\n')
fprintf('nu_res = %.4f THz\n', profiles(1).nu_THz)
fprintf('A_total = %.4f\n', A(k_res))
fprintf('R = %.4f\n', refl(k_res))
fprintf('T = %.4f\n\n', trans(k_res))

for jj = 1:length(layer_names)
    fprintf('%s loss fraction = %.4f\n', layer_names{jj}, loss_frac(jj))
end