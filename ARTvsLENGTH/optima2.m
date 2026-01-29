function refl = optima2(p, thcEXP, nu0) 

    c = 29979245800; % speed of light in cm/s 
    npt = 1;
    nu = nu0;
    omega = 2*pi*nu;  % angular frequency in rad/s
    
    refl = zeros(1, length(thcEXP));

    for jj = 1:length(thcEXP)
        Air = struct('type', 0, ...
                    'd', 1, ...
                    'n0', 1, ...
                    'omega_p', 0, ...
                    'tau', 0, ...
                    'omega_0', 0);

        Au = struct('type', 1, ...
                'd', 2e-4, ...
                'n0', 1, ...
                'omega_p', sqrt(4*pi*9e9*5e7*1e19), ...
                'tau', 1e-19);


        Given = struct('type', -1, ...
                     'd', thcEXP(jj)*1e2, ...
                     'n0', 0, ...
                     'eps1', p(1), ...
                     'eps2', p(2), ...
                     'omega', omega);
        
        structure = {Air, Given, Au, Air};
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
            end
            % n_layers()
            n_layers(j, :) = sqrt(eps);  % complex refractive index
            eps_layers(j, :) = eps;  % complex 
            % n_layers(j, :) = real(n_layers(j, :)) + 1i * abs(g(n_layers(j, :)));
        end
        eps_layers(2) = p(1) + 1j*p(2);
        n_layers(2) = sqrt(eps_layers(2));
        
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
        
        Ein = squeeze(T(:,2,:));  % incoming field from right-side excitation
        
        refl(jj) = abs(Ein(1,:) ./ Ein(2,:)).^2;
        % trans = abs(1 ./ Ein(2,:)).^2;
    end
end
