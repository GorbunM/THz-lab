 
% nthc = 51;
% thcmin = 1e-6; thcmax = 1e-3;
% thclist = exp(linspace(log(thcmin), log(thcmax), nthc));
% npt = 1001;
% numin = 1e11; numax = 2e13;
% nu = linspace(numin, numax, npt);  % frequency in Hz

nthc = 51;
thcmin = 1e-4; thcmax = 3e-3;
thclist = exp(linspace(log(thcmin), log(thcmax), nthc));

npt = 501;
numin = 2e12; numax = 4e13;
nu = linspace(numin, numax, npt);  % frequency in Hz

omega = 2*pi*nu;  % angular frequency in rad/s
Rout = zeros(nthc, npt);
Tout = zeros(nthc, npt);
Aout = zeros(nthc, npt);

%%
for dd = 1:nthc
    
    c = 29979245800; % speed of light in cm/s 
    
    Air = struct('type', 0, ...
                 'd', 1, ...
                 'n0', 1);

    Au = struct('type', 1, ...
                'd', 2e-4, ...
                'n0', 1, ...
                'omega_p', sqrt(4*pi*9e9*5e7*1e19), ...
                'tau', 1e-19);

    MLL = struct('type', 3, ...
                'd', thclist(dd), ...
                'n0', sqrt(3.85), ...
                'omega_p', 4.5e13, ...
                'tauD', 1 / 3.9e12, ...
                'tau', [1 / 3.9e12, 1 / 3e12], ...
                'omega_0', [7.55e13, 6.55e13]);
    
    structure = {Air, MLL, Au, Air};
    N = length(structure);
    
    n_layers = zeros(N, npt);
    eps_layers = zeros(N, npt);
    for j = 1:N
        curlayer = structure{j};
        eps = curlayer.n0 ^ 2;
        if curlayer.type == -1
            eps = curlayer.eps1 + 1j*curlayer.eps2;
        elseif curlayer.type == 1
            w_p = curlayer.omega_p;
            g   = 1 / curlayer.tau;
            eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);              % Drude model
        elseif curlayer.type == 2
            w_p = curlayer.omega_p;
            for peak = 1:length(curlayer.omega_0)
                g   = 1 / curlayer.tau(peak);
                w_0 = curlayer.omega_0(peak);
                eps = eps + (w_p^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Lorentz model
            end
        elseif curlayer.type == 3
            w_p = curlayer.omega_p;
            g   = 1 / curlayer.tauD;
            eps = eps - (w_p^2) ./ (omega.^2 + 1i*g.*omega);
            for peak = 1:length(curlayer.omega_0)
                g   = 1 / curlayer.tau(peak);
                w_0 = curlayer.omega_0(peak);
                eps = eps + (w_p^2) ./ (w_0^2 - omega.^2 - 1i*g.*omega);  % Drude-Lorentz model
            end
            n_layers(j, :) = sqrt(eps);  % complex refractive index
            eps_layers(j, :) = eps;  % complex 
        end
        n_layers(j, :) = sqrt(eps);  % complex refractive index
        eps_layers(j, :) = eps;  % complex refractive index
        % n_layers(j, :) = real(n_layers(j, :)) + 1i * abs(g(n_layers(j, :)));
    end
    
    % Assemble transfer matrices
    
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
    
    % r = (nL - nR) ./ (nL + nR);
    % t = 2 * nL ./ (nL + nR);
    M = zeros(2,2,npt);
    M(1,1,:) = 1 + nL ./ nR;
    M(1,2,:) = 1 - nL ./ nR;
    M(2,1,:) = M(1,2,:);
    M(2,2,:) = M(1,1,:);
    M = 0.5 * M;
    
    T = pagemtimes(T, M);  % final multiply
    
    % 7. Compute reflection and transmission spectra 
    Ein = squeeze(T(:,2,:));  % incoming field from right-side excitation
    
    refl = abs(Ein(1,:) ./ Ein(2,:)).^2;
    trans = abs(1 ./ Ein(2,:)).^2;
    
    A = 1 - refl - trans;  % absorption (if any)
    Rout(dd, :) = refl;
    Tout(dd, :) = trans;
    Aout(dd, :) = A;
end

%%
% figure
% surf(nu, thclist, Tout, 'EdgeColor','none')
% set(gca,'YScale','log');
% xlabel('\nu');
% ylabel('d');
% zlabel('T');

figure
surf(nu, thclist, Rout, 'EdgeColor','none')
set(gca,'YScale','log');
xlabel('\nu');
ylabel('d');
zlabel('R');

% figure
% surf(nu, thclist, Aout, 'EdgeColor','none')
% set(gca,'YScale','log');
% xlabel('\nu');
% ylabel('d');
% zlabel('A');