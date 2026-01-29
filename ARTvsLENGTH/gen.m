%% 1. Constants (CGS units)

c = 29979245800; % speed of light in cm/s 

%% 2. Frequency range

npt = 1;
% nu = linspace(1e11, 2e13, npt);  % frequency in Hz
nu = 2e12;
omega = 2*pi*nu;  % angular frequency in rad/s

%% 3. Layer definitions via Lorentz/Drude parameters

% Layer types: 0 -- air/vacuum
%              1 -- Drude
%              2 -- Lorenz
%              3 -- Drude + Lorenz 

% Layer templates:
    % DrudeLayer = struct('type', 1, ...
    %                     'd', 1e-4, ...
    %                     'n0', 1, ...
    %                     'omega_p', 2*pi*2e14, ...
    %                     'tau', 2*pi*1e13);
    % 
    % LorentzLayer = struct('type', 2, ...
    %                       'd', 5e-5, ...
    %                       'n0', 1, ...
    %                       'omega_p', 2*pi*1e14, ...
    %                       'tau', 2*pi*1e13, ...
    %                       'omega_0', 2*pi*1.5e14);
    %
    % DrudeLorentzLayer = struct('type', 3, ...
    %                       'd', 5e-5, ...
    %                       'n0', 1, ...
    %                       'omega_p', 2*pi*1e14, ...
    %                       'tau', 2*pi*1e13, ...
    %                       'omega_0', 2*pi*1.5e14);
    %
    % DrudeLorentzLayer = struct('type', 31, ...
    %                       'd', 5e-5, ...
    %                       'n0', 1, ...
    %                       'omega_p', 2*pi*1e14, ...
    %                       'tau', 2*pi*1e13, ...
    %                       'omega_0', 2*pi*1.5e14);
% thcEXP = [1e-4, 5e-4, 1e-3, 2e-3, 3e-3];
thcEXP = linspace(1e-6, 3e-5, 5)
refl = zeros(1, length(thcEXP));
for jj = 1:length(thcEXP)
    Air = struct('type', 0, ...
                 'd', 1, ...
                 'n0', 1, ...
                 'omega_p', 0, ...
                 'tau', 0, ...
                 'omega_0', 0);
    
    Si3N4 = struct('type', 2, ...
                'd', thcEXP(jj)*1e2, ...
                'n0', sqrt(3.85), ...
                'omega_p', 4.5e13, ...
                'tau', 1 / 3.9e12, ...
                'omega_0', 7.55e13);
    % disp(Si3N4.d)
    
    Au = struct('type', 1, ...
                'd', 2e-4, ...
                'n0', 1, ...
                'omega_p', sqrt(4*pi*9e9*5e7*1e19), ...
                'tau', 1e-19);
    
    Null = struct('type', 1, ...
                'd', 2e-6, ...
                'n0', 0, ...
                'omega_p', 0, ...
                'tau', 1e-19);
    
    
    
    %Layer 3 gold
    % tauL3=1e-19; thcL3=2*1e-6; sig00L3=1; n03=1;
    
    %% 4. Stack construction
    
    structure = {Air, Si3N4, Au, Air};
    N = length(structure);
    
    %% 5. Compute complex refractive index for each layer at all frequencies
    
    n_layers = zeros(N, npt);
    eps_layers = zeros(N, npt);
    for j = 1:N
        curlayer = structure{j};
        eps = curlayer.n0 ^ 2;
        if curlayer.type == 1
            w_p = curlayer.omega_p;
            g   = 1 / curlayer.tau;
            eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);  % Drude model
        elseif curlayer.type == 2
            w_p = curlayer.omega_p;
            g   = 1 / curlayer.tau;
            w_0 = curlayer.omega_0;
            eps = eps + (w_p^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Lorentz model
        elseif curlayer.type == 3
            w_p = curlayer.omega_p;
            g   = 1 / curlayer.tau;
            w_0 = curlayer.omega_0;
            eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega) + (w_p^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  %  Drude-Lorentz model
        end
        n_layers(j, :) = sqrt(eps);  % complex refractive index
        eps_layers(j, :) = eps;  % complex refractive index
        % n_layers(j, :) = real(n_layers(j, :)) + 1i * abs(g(n_layers(j, :)));
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
    
    r = (nL - nR) ./ (nL + nR);
    t = 2 * nL ./ (nL + nR);
    M = zeros(2,2,npt);
    M(1,1,:) = 1 + nL ./ nR;
    M(1,2,:) = 1 - nL ./ nR;
    M(2,1,:) = M(1,2,:);
    M(2,2,:) = M(1,1,:);
    M = 0.5 * M;
    
    T = pagemtimes(T, M);  % final multiply
    
    %% 7. Compute reflection and transmission spectra 
    Ein = squeeze(T(:,2,:));  % incoming field from right-side excitation
    
    refl(jj) = abs(Ein(1,:) ./ Ein(2,:)).^2;
    trans = abs(1 ./ Ein(2,:)).^2;
    % disp(refl)
    A = 1 - refl - trans;  % absorption (if any)
end
%% Plot

% figure
% hold on
% plot(nu*1e-12, refl, 'r', 'DisplayName','R')
% % plot(nu*1e-12, trans, 'b', 'DisplayName','T')
% xlabel('\nu, THz')
% ylabel('R, T')
% legend
disp(refl)
reflgen = refl;

figure
hold on
scatter(thcEXP, reflgen, 'r', 'DisplayName','R')
% plot(nu*1e-12, trans, 'b', 'DisplayName','T')
% xlabel('\nu, THz')
ylabel('R, T')
legend

