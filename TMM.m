%% 1. Constants (CGS units)

c = 29979245800; % speed of light in cm/s 

%% 2. Frequency range

npt = 500;
nu = linspace(1.5e11, 1.5e13, npt);  % frequency in Hz
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
    %                       'omega_p', 2*pi*1e14, ...
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

Air = struct('type', 0, ...
             'd', 1, ...
             'n0', 1);

Si = struct('type', 1, ...
            'd', 5e-3, ...
            'n0', 3.5, ...
            'omega_p', 0, ...
            'tauD', 1, ...
            'tau', 1, ...
            'omega_0', 1);

Au = struct('type', 1, ...
            'd', 2e-5, ...
            'n0', 1, ...
            'omega_p', sqrt(4*pi*9e9*5e7*1e19), ...
            'tauD', 1e-19);

Lor = struct('type', 2, ...
             'd', 5e-5, ...
             'n0', 1, ...
             'omega_p', 2*pi*1e13, ...
             'tau', 1/(2*pi*1e12), ...
             'omega_0', 2*pi*5e12, ...
             'magn', 1);

LorentzLayer = struct('type', 2, ...
                      'd', 5e-5, ...
                      'n0', 1, ...
                      'omega_p', 2*pi*1e14, ...
                      'tau', [2*pi*1e13, 2*pi*5e12], ...
                      'omega_0', [2*pi*1.5e12, 2*pi*1.2e13], ...
                      'magn', [1, 10]);

%% 4. Stack construction

structure = {Air, LorentzLayer, Air};
N = length(structure);

%% 5. Compute complex refractive index for each layer at all frequencies

n_layers = zeros(N, npt);
eps_layers = zeros(N, npt);
for j = 1:N
    curlayer = structure{j};
    eps = curlayer.n0 ^ 2;
    if curlayer.type == -1
        eps = curlayer.eps1 + 1j*curlayer.eps2;
    elseif curlayer.type == 1
        w_p = curlayer.omega_p;
        g   = 1 / curlayer.tauD;
        eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);              % Drude model
    elseif curlayer.type == 2
        w_p = curlayer.omega_p;
        for peak = 1:length(curlayer.omega_0)
            g   = 1 / curlayer.tau(peak);
            w_0 = curlayer.omega_0(peak);
            A_0 = curlayer.magn(peak);
            eps = eps + A_0 * (w_p^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Lorentz model
        end
    elseif curlayer.type == 3
        w_p = curlayer.omega_p;
        A_0 = curlayer.magn;
        g   = 1 / curlayer.tauD;
        eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);
        for peak = 1:length(curlayer.omega_0)
            g   = 1 / curlayer.tau(peak);
            w_0 = curlayer.omega_0(peak);
            eps = eps + A_0(peak) * (w_p^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Drude-Lorentz model
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

A = 1 - refl - trans;  % absorption (if any)

%% Plot

nu = nu*1e-12;

figure
hold on
% plot(nu, refl, 'r', 'DisplayName','R')
% plot(nu, trans, 'b', 'DisplayName','T')
plot(nu, real(eps_layers(2,:)), 'DisplayName', '$\varepsilon^{\prime}$')
plot(nu, imag(eps_layers(2,:)), 'DisplayName', '$\varepsilon^{\prime\prime}$')
xlabel('\nu, THz')
% ylabel('R, T')
legend('Interpreter','latex', 'FontSize', 14)