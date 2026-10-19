%% 1. Constants (CGS units)

c = 29979245800; % speed of light in cm/s 

%% 2. Frequency range

npt = 2000;
nu = linspace(1.5e11, 30e12, npt);  % frequency in Hz

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
             'd', 2e-5, ...
             'n0', 1);

    Di = struct('type', 1, ...
                'd', 5e-2, ...
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
                'omega_0', 2*pi*10.4e12, ... %2*pi*1.05e13, ...
                'magn', 4);

    CNT = struct('type', 1, ...
                'd', 2e-6, ...
                'n0', 1, ...
                'omega_p', Omp, ...
                'tauD', 1e-13);
    
    CNT1 = struct('type', 1, ...
                'd', 6*2e-6, ...
                'n0', 1, ...
                'omega_p', Omp, ...
                'tauD', 1e-13);

%% 4. Stack construction

% structure = {Air, CNT1, Air};
structure = {Air, CNT, Air, CNT, Air, CNT, Air, CNT, Air, CNT, Air, CNT, Air};
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
Dd = [];
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
    Dd(j, :) = delta;

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
complexT = 1 ./ Ein(2,:);

A = 1 - refl - trans;  % absorption (if any)

%% Plot

nuu = nu*1e-12;

figure
hold on
plot(nuu, A, 'k', 'DisplayName','A', 'LineWidth', 2)
plot(nuu, refl, 'r', 'DisplayName','R','LineWidth', 2)
plot(nuu, trans, 'b', 'DisplayName','T','LineWidth', 2)
ylim([0, 1])

% plot(nuu, phase, 'DisplayName','$\phi$')
% ylim([0, pi/2])

% plot(nuu, real(eps_layers(2,:)), 'DisplayName', '$\varepsilon^{\prime}$')
% plot(nuu, imag(eps_layers(2,:)), 'DisplayName', '$\varepsilon^{\prime\prime}$')

% xlabel('\nu, THz')
% ylabel('Emissivity')
ax = gca;
ax.FontSize = 14;

legend('Interpreter','latex', 'FontSize', 24, 'Location','northeast')

box on
grid on

% %% 8. Field and loss profiles from TMM
% 
% % Frequencies for comparison
% [~, k_res] = max(A);              % resonance frequency
% [~, k_off] = min(abs(nu - 8e12)); % off-resonance reference, change if needed
% 
% k_list = [k_res, k_off];
% labels = {'resonance', 'off-resonance'};
% 
% profiles = struct([]);
% 
% for kk = 1:length(k_list)
% 
%     k_ind = k_list(kk);
% 
%     % Complex reflection amplitude, not reflectance
%     r = Ein(1,k_ind) / Ein(2,k_ind);
% 
%     % Initial amplitudes in incident medium:
%     % E(z) = E_plus exp(i beta z) + E_minus exp(-i beta z)
%     amp = [1; r];
% 
%     z_all = [];
%     E2_all = [];
%     ploss_all = [];
%     layer_id = [];
% 
%     z_offset = 0;
% 
%     for j = 2:N-1
% 
%         % Interface from layer j-1 to layer j
%         nL = n_layers(j-1,k_ind);
%         nR = n_layers(j,k_ind);
% 
%         M = 0.5 * [1 + nL/nR, 1 - nL/nR; ...
%                    1 - nL/nR, 1 + nL/nR];
% 
%         % Amplitudes just inside layer j
%         amp = M \ amp;
% 
%         % Local coordinate inside layer j
%         d = structure{j}.d;
% 
%         dz = 1e-7;                 % spatial step in cm, 1 nm
%         Nz = ceil(d/dz);  % number of points in current layer
% 
%         zloc = linspace(0, d, Nz);
% 
%         beta = omega(k_ind) * nR / c;
% 
%         E = amp(1) * exp(1i * beta * zloc) + ...
%             amp(2) * exp(-1i * beta * zloc);
% 
%         E2 = abs(E).^2;
% 
%         % CGS loss density up to physical normalization:
%         % q(z) = omega/(8*pi) * Im(eps) * |E|^2
%         eps_im = imag(eps_layers(j,k_ind));
%         ploss = omega(k_ind) / (8*pi) * eps_im * E2;
% 
%         z_all = [z_all, z_offset + zloc];
%         E2_all = [E2_all, E2];
%         ploss_all = [ploss_all, ploss];
%         layer_id = [layer_id, j * ones(size(zloc))];
% 
%         % Propagate amplitudes to the right boundary of layer j
%         P = [exp(1i * beta * d), 0; ...
%              0, exp(-1i * beta * d)];
% 
%         amp = P * amp;
% 
%         z_offset = z_offset + d;
%     end
% 
%     profiles(kk).k_ind = k_ind;
%     profiles(kk).nu_THz = nu(k_ind) * 1e-12;
%     profiles(kk).z = z_all;
%     profiles(kk).ploss = ploss_all;
%     profiles(kk).layer_id = layer_id;
% end
% 
% 
% %% 11. Power-loss density profile at resonance
% 
% figure
% plot(profiles(1).z * 1e4, profiles(1).ploss, 'k', 'LineWidth', 2)
% hold on
% 
% xline(z_TL, ':', 'TL/Di', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom')
% xline(z_Di, ':', 'Di/Au', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom')
% 
% xlabel('z, \mum')
% ylabel('q(z), arb. units')
% title(['Power-loss density at \nu = ', num2str(profiles(1).nu_THz, '%.2f'), ' THz'])
% 
% ax = gca;
% ax.FontSize = 24;
% box on
% grid on
%%
% Phase difference

% Rr = (n_layers(1:4) - n_layers(2:5))./(n_layers(1:4) + n_layers(2:5));
% Tf = (2*n_layers(1:4))./(n_layers(1:4) + n_layers(2:5));
% Tb = (2*n_layers(2:5))./(n_layers(1:4) + n_layers(2:5));
% 
% rr = unwrap(angle(Rr));
% tf = unwrap(angle(Tf));
% tb = unwrap(angle(Tb));
% dd = unwrap(angle(Dd));
% 
% ksi = tf'+tb'+2*dd;
% 
% phi0 = rr(1);
% phi1 = ksi(1, :) + rr(2);
% phi2 = ksi(1, :) + ksi(2, :) + rr(3);
% phi3 = ksi(1, :) + ksi(2, :) + ksi(3, :) + rr(4);
% 
% Phi = [phi1; phi2; phi3] - phi0;
% Phi = Phi/pi;
% figure
% hold on
% for k = 1:3
%     plot(nuu, Phi(k, :));
% end
% 
% legend