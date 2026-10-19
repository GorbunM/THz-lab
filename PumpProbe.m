%% 0. Settings and data loading

clear; clc; 

curdir = "C:\Users\mgorbun\Documents\MATLAB\pump-probe\2026-10-15 670 nm\parralel";  % Data directory
frmin = 0.15;
frmax = 1.8;
% frmax = -1;

files = dir(fullfile(curdir, '*.dat'));

TT = {};
names = {};
delay = [];
type = {};

for i = 1:numel(files)

    name = files(i).name;
    lowname = lower(name);

    if startsWith(lowname, 'reference')
        kind = 'reference';
        dt = NaN;

    elseif startsWith(lowname, 'nopump')
        kind = 'nopump';
        dt = NaN;

    elseif startsWith(lowname, 'pls') || startsWith(lowname, 'min')
        kind = 'pump';
    
        token = regexp(lowname, ...
            '^(pls|min)(\d+(?:\.\d+)?)ps\d*\.dat$', ...
            'tokens', 'once');
    
        if isempty(token), continue; end
    
        str = token{2};
        dt = str2double(str);
    
        % 025 -> 0.25, 05 -> 0.5, 075 -> 0.75
        if startsWith(str, '0') && numel(str) > 1 && ~contains(str, '.')
            dt = dt / 10^(numel(str)-1);
        end
    
        if strcmp(token{1}, 'min')
            dt = -dt;
        end

    else
        continue;
    end

    data = readmatrix(fullfile(curdir, name));

    TT{end+1} = data(:,1:2);
    names{end+1} = name;
    delay(end+1) = dt;
    type{end+1} = kind;
end

% Indices
idxRef = find(strcmp(type, 'reference'));
idxNoPump = find(strcmp(type, 'nopump'));
idxPump = find(strcmp(type, 'pump'));

% Sort pump measurements by delay
[delayPump, order] = sort(delay(idxPump));
idxPump = idxPump(order);

% Fourier transforms
FFT = cell(size(TT));
S = cell(size(TT));

for i = 1:numel(TT)

    FFT{i} = TDS_FFT(TT{i}, frmin, frmax);

    S{i} = FFT{i}(:,2) .* exp(1i*FFT{i}(:,3));
end

freq = FFT{1}(:,1);

% Average complex spectra
Sref = mean(cat(2, S{idxRef}), 2);
SnoPump = mean(cat(2, S{idxNoPump}), 2);

% Average time traces
TTref = mean(cat(3, TT{idxRef}), 3);
TTnoPump = mean(cat(3, TT{idxNoPump}), 3);

% Normalization
Href = cell(size(S));
HnoPump = cell(size(S));

for i = 1:numel(S)
    Href{i} = S{i} ./ Sref;
    HnoPump{i} = S{i} ./ SnoPump;
end

% Colors
colors = turbo(max(numel(idxPump), 2));

%% 1-8. All plots in subplots

figure('Name','TDS overview');
set(gcf, 'Position', [100 100 1800 900]);
tiledlayout('flow');
%% 1. All time traces
nexttile;
hold on;

% References (dashed)
for j = 1:numel(idxRef)
    i = idxRef(j);
    plot(TT{i}(:,1), TT{i}(:,2), '--', ...
        'DisplayName', names{i});
end

% No pump
for j = 1:numel(idxNoPump)
    i = idxNoPump(j);
    plot(TT{i}(:,1), TT{i}(:,2), ...
        'DisplayName', names{i});
end

% Pumped
for j = 1:numel(idxPump)
    i = idxPump(j);
    plot(TT{i}(:,1), TT{i}(:,2), ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('%g ps', delay(i)));
end

xlabel('Time (ps)');
ylabel('Signal (a.u.)');
title('Time traces');
grid on;
legend('Interpreter','none','Location','eastoutside');


%% 2. All Fourier spectra
nexttile;
hold on;

% References (dashed)
for j = 1:numel(idxRef)
    i = idxRef(j);
    plot(freq, abs(S{i}), '--', ...
        'DisplayName', names{i});
end

% No pump
for j = 1:numel(idxNoPump)
    i = idxNoPump(j);
    plot(freq, abs(S{i}), ...
        'DisplayName', names{i});
end

% Pumped
for j = 1:numel(idxPump)
    i = idxPump(j);
    plot(freq, abs(S{i}), ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('%g ps', delay(i)));
end

xlabel('Frequency (THz)');
ylabel('|FFT| (a.u.)');
title('Fourier spectra');
grid on;
legend('Interpreter','none','Location','eastoutside');


%% 3. Amplitude normalized to reference
nexttile;
hold on;

plot(freq, abs(SnoPump ./ Sref), ...
    'k', 'LineWidth', 2, ...
    'DisplayName', 'No pump');

for j = 1:numel(idxPump)
    i = idxPump(j);
    plot(freq, abs(Href{i}), ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('%g ps', delay(i)));
end

xlabel('Frequency (THz)');
ylabel('|E / E_{ref}|');
title('Amp / ref');
grid on;
legend('Location','eastoutside');


%% 4. Phase normalized to reference
nexttile;
hold on;

plot(freq, unwrap(angle(SnoPump ./ Sref)), ...
    'k', 'LineWidth', 2, ...
    'DisplayName', 'No pump');

for j = 1:numel(idxPump)
    i = idxPump(j);
    plot(freq, unwrap(angle(Href{i})), ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('%g ps', delay(i)));
end

xlabel('Frequency (THz)');
ylabel('Phase (rad)');
title('Phase / ref');
grid on;
legend('Location','eastoutside');


%% 5. Amplitude normalized to no pump
% nexttile;
% hold on;
% 
% yline(1, 'k--', 'No pump');
% 
% for j = 1:numel(idxPump)
%     i = idxPump(j);
%     plot(freq, abs(HnoPump{i}), ...
%         'Color', colors(j,:), ...
%         'DisplayName', sprintf('%g ps', delay(i)));
% end
% 
% xlabel('Frequency (THz)');
% ylabel('|E / E_{no pump}|');
% title('Amp / no pump');
% grid on;
% legend('Location','eastoutside');
% 
% 
%% 6. Phase normalized to no pump
% nexttile;
% hold on;
% 
% yline(0, 'k--', 'No pump');
% 
% for j = 1:numel(idxPump)
%     i = idxPump(j);
%     plot(freq, unwrap(angle(HnoPump{i})), ...
%         'Color', colors(j,:), ...
%         'DisplayName', sprintf('%g ps', delay(i)));
% end
% 
% xlabel('Frequency (THz)');
% ylabel('Phase (rad)');
% title('Phase / no pump');
% grid on;
% legend('Location','eastoutside');


%% 7. Peak-to-peak amplitude vs delay
Amp = zeros(size(idxPump));

for j = 1:numel(idxPump)
    i = idxPump(j);
    signal = TT{i}(:,2);
    Amp(j) = max(signal) - min(signal);
end

Amp0 = zeros(size(idxNoPump));

for j = 1:numel(idxNoPump)
    signal = TT{idxNoPump(j)}(:,2);
    Amp0(j) = max(signal) - min(signal);
end

nexttile;
hold on;
plot(delayPump, Amp, 'o-', 'LineWidth', 1.5, 'DisplayName','Pump');
yline(mean(Amp0), 'k--', 'No pump');
xlabel('Pump-probe delay (ps)');
ylabel('Peak-to-peak amp (a.u.)');
title('TT amplitude vs delay');
grid on;


%% 8. Integrated squared signal vs delay
Energy = zeros(size(idxPump));

for j = 1:numel(idxPump)
    i = idxPump(j);
    t = TT{i}(:,1);
    signal = TT{i}(:,2);
    Energy(j) = trapz(t, abs(signal).^2);
end

Energy0 = zeros(size(idxNoPump));

for j = 1:numel(idxNoPump)
    i = idxNoPump(j);
    Energy0(j) = trapz(TT{i}(:,1), abs(TT{i}(:,2)).^2);
end

nexttile;
hold on;
plot(delayPump, Energy, 'o-', 'LineWidth', 1.5, 'DisplayName','Pump');
yline(mean(Energy0), 'k--', 'No pump');
xlabel('Pump-probe delay (ps)');
ylabel('\int |E|^2 dt (a.u.^2 ps)');
title('Integrated TT intensity');
grid on;

sgtitle('TDS data overview');