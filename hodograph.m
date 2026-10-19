%% 1. Constants (CGS units)

c = 29979245800; % speed of light in cm/s 

% A0 = readmatrix("C:\Users\mgorbun\Documents\MATLAB\TDS\setup-test\2025-0224-empty-LG.csv");
A0 = readmatrix("C:\Users\mgorbun\Documents\MATLAB\sample10-2\air_0deg.txt");

x = A0(:,1);
y = A0(:,3);

xNew = linspace(x(1), x(end), 1*numel(x)).';
yNew = interp1(x, y, xNew, 'pchip');


dt = mean(diff(xNew));
Nfft = 8*numel(yNew);  % нечётная длина, окно около 160 пс

tFFT = xNew(1) + (0:Nfft-1).' * dt;
yFFT = zeros(Nfft,1);
yFFT(1:numel(yNew)) = yNew;

TT = [tFFT, yFFT];
% TTdiffnorm = [xNew(1:end-1), diff(normpdf(xNew, xNew(floor(length(xNew)/15)), 0.05))];
% TT = TTdiffnorm;

Ref = traceFFT(TT);
%% 2. Frequency range

npt = length(Ref);
nu = Ref(:,1).' * 1e12;

Rsh = 400; % Ohm/sq
d = 2e-6; %cm
tau = 1e-13; %s

Omp = 10*sqrt((4*pi*9e9)/(Rsh*d*tau));

omega = 2*pi*nu;  % angular frequency in rad/s


    Air = struct('type', 0, ...
             'd', 1, ...
             'n0', 1);

    CNT = struct('type', 1, ...
                'd', 2e-6, ...
                'n0', 1, ...
                'omega_p', Omp, ...
                'tauD', 1e-13);
    
    Di = struct('type', 0, ...
                'd', 5e-2, ...
                'n0', 3.5);

%% 4. Stack construction

structure = {Air, CNT, Di};
% structure = {Air, Air};
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
Tpar = 1 ./ Ein(2,:);


%% 4. Stack construction

structure = {Air, Di};
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
Tper = 1 ./ Ein(2,:);

dc = (omega == 0);
Tpar(dc) = 0;
Tper(dc) = 0;

Epar = Ref(:,2).*Tpar';
Eper = Ref(:,2).*Tper';

Epart = real(ifft(Epar));
Epert = real(ifft(Eper));

magn = max(max(Epert),max(Epart));

figure
tiledlayout(3, 2)

plotidx = (1:floor((length(xNew)-1)*1));

nexttile([3 1])
plot(Epart(plotidx), Epert(plotidx))
title('Hodograph')
axis equal

% Справа сверху
nexttile
plot(xNew(plotidx), Epart(plotidx))
title('parallel')
% xlim([0, numel(yNew)])

% Справа снизу
nexttile
plot(xNew(plotidx), Epert(plotidx))
title('perpendicular')
% xlim([0, numel(yNew)])

nexttile
% plot(yNew)
plot(xNew(plotidx), TT(plotidx,2))
title('initial')


%%
% v = VideoWriter('R.mp4','MPEG-4');
% v.FrameRate = 5;
% open(v)
% figure
% % for k = 1:125
% for k = 1:numel(Epart)
%     xxx = Epart(1:k);
%     yyy = Epert(1:k);
% 
%     tiledlayout(3, 2)
% 
%     nexttile([3 1])
%     plot(xxx, yyy)
%     title('Hodograph')
%     axis equal
% 
%     nexttile
%     plot(xxx)
%     title('parallel')
% 
%     nexttile
%     plot(yyy)
%     title('perpendicular')
% 
%     nexttile
%     % plot(yNew)
%     plot(TT(1:k,2))
%     title('initial')
% 
%     drawnow
%     writeVideo(v, getframe(gcf));
% end
% 
% curdir = "C:\Users\mgorbun\Documents\MATLAB\THz-lab";
% saveas(gcf, fullfile(curdir, 'Japan-interp-maxf.png'))
% 
% close(v)

%%

% el = abs(abs(Tpar) - abs(Tper))./(abs(Tpar) + abs(Tper));
% plot(el)

%%

% figure
% title('Hodograph')
% 
% hold on
% out = [];
% for theta = (0:4)*pi/8
%     plot(Epart(plotidx)./magn*cos(theta), Epert(plotidx)./magn*sin(theta))
% end
% axis equal

%%
% figure
% xlabel('THz')
% ylabel('changeel')
% grid on
% hold on
% 
% for theta = (0:4)*pi/8
% 
%     TTelpar = [xNew(plotidx), Epart(plotidx)*cos(theta)];
%     TTelper = [xNew(plotidx), Epert(plotidx)*sin(theta)];
% 
% 
%     Tel0 = traceFFT([xNew(plotidx), yNew(plotidx)]);
%     Telpar = traceFFT(TTelpar);
%     Telper = traceFFT(TTelper);
% 
%     Ecra0 = abs((Tel0(:,2) + 1j * Tel0(:,2)) / sqrt(2));
%     Ecri0 = abs((Tel0(:,2) - 1j * Tel0(:,2)) / sqrt(2));
% 
%     Ecra = abs((Telpar(:,2) + 1j * Telper(:,2)) / sqrt(2));
%     Ecri = abs((Telpar(:,2) - 1j * Telper(:,2)) / sqrt(2));
% 
%     changeel = asin(abs(Ecri - Ecra)./(Ecra + Ecri))./pi*180;
%     changeel0 = asin(abs(Ecri0 - Ecra0)./(Ecra0 + Ecri0))./pi*180;
% 
%     idx = (0 < Telper(:,1)) & (Telper(:,1) < 2.5);
% 
%     plot(Telper(idx,1), changeel(idx,1)-changeel0(idx,1), 'LineWidth', 3)
% 
% end
% legend('FontSize', 24)
% disp('That''s it')

%%

% folder = "C:\Users\mgorbun\Documents\MATLAB\hodograph\";
% 
% hold on
% for theta = (0:4)*pi/8
%     Ex = Epart(plotidx)./magn*cos(theta);
%     Ey = Epert(plotidx)./magn*sin(theta);
% 
%     plot(Ex, Ey)
% 
%     filename = sprintf('theta_%g_deg.csv', rad2deg(theta));
%     writematrix([Ex(:), Ey(:)], fullfile(folder, filename));
% end
% axis equal

%%

% folder = "C:\Users\mgorbun\Documents\MATLAB";
% 
% figure
% xlabel('THz')
% ylabel('changeel')
% grid on
% hold on
% 
% Tel0 = traceFFT([xNew(plotidx), yNew(plotidx)]);
% 
% Ecra0 = abs((Tel0(:,2) + 1j * Tel0(:,2)) / sqrt(2));
% Ecri0 = abs((Tel0(:,2) - 1j * Tel0(:,2)) / sqrt(2));
% 
% changeel0 = asin(abs(Ecri0 - Ecra0)./(Ecra0 + Ecri0))/pi*180;
% 
% fflag = 0;
% linest = {'-', ':', '--'};
% for theta = (1:7)*pi/16
%     fflag = fflag+1
%     flag = floor(fflag/5.5)+1
%     ls = linest{flag};
%     TTelpar = [xNew(plotidx), Epart(plotidx)*cos(theta)];
%     TTelper = [xNew(plotidx), Epert(plotidx)*sin(theta)];
% 
%     Telpar = traceFFT(TTelpar);
%     Telper = traceFFT(TTelper);
% 
%     Ecra = abs((Telpar(:,2) + 1j * Telper(:,2)) / sqrt(2));
%     Ecri = abs((Telpar(:,2) - 1j * Telper(:,2)) / sqrt(2));
% 
%     % changeel = asin(abs(Ecri - Ecra)./(Ecra + Ecri))/pi*180;
%     changeel = abs(Ecri - Ecra)./(Ecra + Ecri);
% 
%     idx = (0 < Telper(:,1)) & (Telper(:,1) < 2.5);
% 
%     frequency = Telper(idx,1);
%     delta = changeel(idx,1) - changeel0(idx,1);
% 
%     plot(frequency, delta, ls, 'LineWidth', 3, ...
%         'DisplayName', sprintf('%g°', rad2deg(theta)));
% 
%     filename = sprintf('changeel_theta_%g_deg.csv', rad2deg(theta));
%     writematrix([frequency, delta], fullfile(folder, filename));
% end
% 
% legend('FontSize', 24)
% disp('That''s it')


%%
figure
xlabel('THz')
ylabel('changeel')
grid on
hold on

for theta = (3:3)*pi/8

    Telpar = Tpar*cos(theta);
    Telper = Tper*sin(theta);
   
    x0 = abs(Telpar);
    y0 = abs(Telper);

    Ecra = abs((Telpar + 1j * Telper) / sqrt(2));
    Ecri = abs((Telpar - 1j * Telper) / sqrt(2));
    
    changeel = asin(abs(Ecri - Ecra)./(Ecra + Ecri))./pi*180;
    deltaPhi = angle(Telpar./Telper);

    Aa = x0.^2 + y0.^2;
    
    Rr = x0.^2 - y0.^2;
    Kk = 2*cos(deltaPhi).*x0.*y0;
    Sq = sqrt(Rr.^2 + Kk.^2);
    
    el = sqrt((Aa - Sq)./(Aa + Sq));
    el_angle = asind(el);

    idx = (nu > 0) & (nu < 2.5e12) &  (mod(reshape(1:numel(nu), size(nu)), 10) == 0);
    plot(nu(idx), changeel(idx), 'LineWidth', 3)
    % plot(nu(idx), el_angle(idx), ':', 'LineWidth', 3)
    % plot(nu(idx), deltaPhi(idx), 'LineWidth', 3)

end
legend('FontSize', 24)
disp('That''s it')
