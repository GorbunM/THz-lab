%% 0. Settings and data loading

clear; clc; close all;

flag_save_csv = 1;
flag_save_png = 1;

curdir = "C:\Users\mgorbun\Documents\MATLAB\2026-10-30 cnt90 1113 nm repeat with less gap\2026-10-30 cnt90 1113 nm repeat with less gap\perpendicular 1113 nm";  % Data directory
frmin = 0.15;
frmax = 1.3;
% frmax = -1;

saveDir = fullfile(curdir, 'export');
pngResolution = 300;

if flag_save_csv || flag_save_png
    if ~isfolder(saveDir)
        mkdir(saveDir);
    end
end

[~, dataTag] = fileparts(char(curdir));
dataTag = regexprep(dataTag, '[^a-zA-Z0-9_-]', '_');
overviewTag = sprintf('%s_frmin_%g_frmax_%gTHz', dataTag, frmin, frmax);

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

figOverview = figure('Name','TDS overview');
set(gcf, 'Position', [100 100 1800 900]);
tiledlayout('flow');
% %% 1. All time traces
% nexttile;
% hold on;
% 
% % References (dashed)
% for j = 1:numel(idxRef)
%     i = idxRef(j);
%     plot(TT{i}(:,1), TT{i}(:,2), '--', ...
%         'DisplayName', names{i});
% end
% 
% % No pump
% for j = 1:numel(idxNoPump)
%     i = idxNoPump(j);
%     plot(TT{i}(:,1), TT{i}(:,2), ...
%         'DisplayName', names{i});
% end
% 
% % Pumped
% for j = 1:numel(idxPump)
%     i = idxPump(j);
%     plot(TT{i}(:,1), TT{i}(:,2), ...
%         'Color', colors(j,:), ...
%         'DisplayName', sprintf('%g ps', delay(i)));
% end
% 
% xlabel('Time (ps)');
% ylabel('Signal (a.u.)');
% title('Time traces');
% grid on;
% legend('Interpreter','none','Location','eastoutside');
% 
% 
% %% 2. All Fourier spectra
% nexttile;
% hold on;
% 
% % References (dashed)
% for j = 1:numel(idxRef)
%     i = idxRef(j);
%     plot(freq, abs(S{i}), '--', ...
%         'DisplayName', names{i});
% end
% 
% % No pump
% for j = 1:numel(idxNoPump)
%     i = idxNoPump(j);
%     plot(freq, abs(S{i}), ...
%         'DisplayName', names{i});
% end
% 
% % Pumped
% for j = 1:numel(idxPump)
%     i = idxPump(j);
%     plot(freq, abs(S{i}), ...
%         'Color', colors(j,:), ...
%         'DisplayName', sprintf('%g ps', delay(i)));
% end
% 
% xlabel('Frequency (THz)');
% ylabel('|FFT| (a.u.)');
% title('Fourier spectra');
% grid on;
% legend('Interpreter','none','Location','eastoutside');


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


% %% 5. Amplitude normalized to no pump
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
% %% 6. Phase normalized to no pump
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
% 
% 
% %% 7. Peak-to-peak amplitude vs delay
% Amp = zeros(size(idxPump));
% 
% for j = 1:numel(idxPump)
%     i = idxPump(j);
%     signal = TT{i}(:,2);
%     Amp(j) = max(signal) - min(signal);
% end
% 
% Amp0 = zeros(size(idxNoPump));
% 
% for j = 1:numel(idxNoPump)
%     signal = TT{idxNoPump(j)}(:,2);
%     Amp0(j) = max(signal) - min(signal);
% end
% 
% nexttile;
% hold on;
% plot(delayPump, Amp, 'o-', 'LineWidth', 1.5, 'DisplayName','Pump');
% yline(mean(Amp0), 'k--', 'No pump');
% xlabel('Pump-probe delay (ps)');
% ylabel('Peak-to-peak amp (a.u.)');
% title('TT amplitude vs delay');
% grid on;
% 
% 
% %% 8. Integrated squared signal vs delay
% Energy = zeros(size(idxPump));
% 
% for j = 1:numel(idxPump)
%     i = idxPump(j);
%     t = TT{i}(:,1);
%     signal = TT{i}(:,2);
%     Energy(j) = trapz(t, abs(signal).^2);
% end
% 
% Energy0 = zeros(size(idxNoPump));
% 
% for j = 1:numel(idxNoPump)
%     i = idxNoPump(j);
%     Energy0(j) = trapz(TT{i}(:,1), abs(TT{i}(:,2)).^2);
% end
% 
% nexttile;
% hold on;
% plot(delayPump, Energy, 'o-', 'LineWidth', 1.5, 'DisplayName','Pump');
% yline(mean(Energy0), 'k--', 'No pump');
% xlabel('Pump-probe delay (ps)');
% ylabel('\int |E|^2 dt (a.u.^2 ps)');
% title('Integrated TT intensity');
% grid on;
% 
% sgtitle('TDS data overview');

%% Save the eight overview datasets and one overview image

if flag_save_csv
    % 1. Each trace has its own time column; unequal lengths are padded with NaN.
    columns = {};
    headers = {};
    plotIndices = [idxRef(:); idxNoPump(:); idxPump(:)].';
    for i = plotIndices
        [~, traceName] = fileparts(names{i});
        traceTag = sprintf('trace%d_%s', i, traceName);
        columns(end+1:end+2) = {TT{i}(:,1), TT{i}(:,2)};
        headers(end+1:end+2) = {[traceTag '_time_ps'], [traceTag '_signal_au']};
    end
    saveColumnsCsv(fullfile(saveDir, ['01_time_traces_' overviewTag '.csv']), ...
        columns, headers);

    % 2. Fourier amplitudes.
    columns = {freq};
    headers = {'frequency_THz'};
    for i = plotIndices
        [~, traceName] = fileparts(names{i});
        columns{end+1} = abs(S{i});
        headers{end+1} = sprintf('trace%d_%s_FFT_amplitude_au', i, traceName);
    end
    saveColumnsCsv(fullfile(saveDir, ['02_fourier_spectra_' overviewTag '.csv']), ...
        columns, headers);

    % 3-6. Transmission curves, including the no-pump baseline.
    ampRefColumns = {freq, abs(SnoPump ./ Sref)};
    phaseRefColumns = {freq, unwrap(angle(SnoPump ./ Sref))};
    ampNoPumpColumns = {freq, ones(size(freq))};
    phaseNoPumpColumns = {freq, zeros(size(freq))};
    ampHeaders = {'frequency_THz', 'no_pump_amplitude_ratio'};
    phaseHeaders = {'frequency_THz', 'no_pump_phase_rad'};

    for j = 1:numel(idxPump)
        i = idxPump(j);
        traceTag = sprintf('pump%d_delay_%gps', j, delayPump(j));
        ampRefColumns{end+1} = abs(Href{i});
        phaseRefColumns{end+1} = unwrap(angle(Href{i}));
        ampNoPumpColumns{end+1} = abs(HnoPump{i});
        phaseNoPumpColumns{end+1} = unwrap(angle(HnoPump{i}));
        ampHeaders{end+1} = [traceTag '_amplitude_ratio'];
        phaseHeaders{end+1} = [traceTag '_phase_rad'];
    end

    saveColumnsCsv(fullfile(saveDir, ['03_amplitude_ref_' overviewTag '.csv']), ...
        ampRefColumns, ampHeaders);
    saveColumnsCsv(fullfile(saveDir, ['04_phase_ref_' overviewTag '.csv']), ...
        phaseRefColumns, phaseHeaders);
    saveColumnsCsv(fullfile(saveDir, ['05_amplitude_no_pump_' overviewTag '.csv']), ...
        ampNoPumpColumns, ampHeaders);
    saveColumnsCsv(fullfile(saveDir, ['06_phase_no_pump_' overviewTag '.csv']), ...
        phaseNoPumpColumns, phaseHeaders);

    % 7-8. Delay dependence and the plotted baseline.
    saveColumnsCsv(fullfile(saveDir, ['07_peak_to_peak_' overviewTag '.csv']), ...
        {delayPump, Amp, repmat(mean(Amp0), size(delayPump))}, ...
        {'delay_ps', 'pump_peak_to_peak_au', 'no_pump_peak_to_peak_au'});
    saveColumnsCsv(fullfile(saveDir, ['08_integrated_intensity_' overviewTag '.csv']), ...
        {delayPump, Energy, repmat(mean(Energy0), size(delayPump))}, ...
        {'delay_ps', 'pump_integrated_intensity_au2_ps', ...
        'no_pump_integrated_intensity_au2_ps'});
end

if flag_save_png
    exportgraphics(figOverview, ...
        fullfile(saveDir, ['overview_' overviewTag '.png']), ...
        'Resolution', pngResolution);
end





% %% 9. Interactive transmission vs delay
% 
% % Prepare transmission matrices
% Tref = cat(2, Href{idxPump});
% TnoPump = cat(2, HnoPump{idxPump});
% 
% Tbaseline = SnoPump ./ Sref;
% 
% % Figure
% fig = uifigure('Name', 'Transmission vs delay', ...
%     'Position', [100 100 1200 700]);
% 
% axSpec = uiaxes(fig, 'Position', [70 310 500 340]);
% axDelay = uiaxes(fig, 'Position', [650 310 500 340]);
% 
% % Frequency slider
% uilabel(fig, 'Position', [70 255 250 22], ...
%     'Text', 'Central frequency (THz)');
% 
% fmin = min(freq);
% fmax = max(freq);
% 
% fTicks = ceil(fmin/0.5)*0.5 : 0.5 : floor(fmax/0.5)*0.5;
% fMinor = ceil(fmin/0.1)*0.1 : 0.1 : floor(fmax/0.1)*0.1;
% 
% fSlider = uislider(fig, ...
%     'Position', [70 235 470 3], ...
%     'Limits', [fmin fmax], ...
%     'Value', min(max(1, fmin), fmax), ...
%     'MajorTicks', fTicks, ...
%     'MinorTicks', fMinor);
% 
% % Bandwidth slider
% uilabel(fig, 'Position', [70 175 250 22], ...
%     'Text', 'Averaging bandwidth (THz)');
% 
% bwMax = fmax - fmin;
% 
% bwTicks = 0:0.5:bwMax;
% bwMinor = 0:0.1:bwMax;
% 
% bwSlider = uislider(fig, ...
%     'Position', [70 155 470 3], ...
%     'Limits', [0 bwMax], ...
%     'Value', min(0.1, bwMax), ...
%     'MajorTicks', bwTicks, ...
%     'MinorTicks', bwMinor);
% 
% % Normalization
% uilabel(fig, 'Position', [650 255 200 22], ...
%     'Text', 'Normalization');
% 
% normMenu = uidropdown(fig, ...
%     'Position', [650 225 210 30], ...
%     'Items', {'Reference', 'No pump'});
% 
% % Component
% uilabel(fig, 'Position', [900 255 200 22], ...
%     'Text', 'Component');
% 
% compMenu = uidropdown(fig, ...
%     'Position', [900 225 210 30], ...
%     'Items', {'Amplitude', 'Real', 'Imaginary', 'Phase'});
% 
% % Value labels
% fLabel = uilabel(fig, ...
%     'Position', [70 275 250 20]);
% 
% bwLabel = uilabel(fig, ...
%     'Position', [70 195 250 20]);
% 
% % Fit result
% tauLabel = uilabel(fig, ...
%     'Position', [650 155 480 35], ...
%     'Text', '', ...
%     'FontSize', 16, ...
%     'FontWeight', 'bold');
% 
% % Update callback
% refresh = @(f0,bw) updateTransmissionDelay( ...
%     axSpec, axDelay, freq, delayPump, ...
%     Tref, TnoPump, Tbaseline, ...
%     f0, bw, normMenu.Value, compMenu.Value, ...
%     fLabel, bwLabel, tauLabel);
% 
% % Continuous updates while dragging
% fSlider.ValueChangingFcn = @(~,event) ...
%     refresh(round(event.Value,2), round(bwSlider.Value,2));
% 
% bwSlider.ValueChangingFcn = @(~,event) ...
%     refresh(round(fSlider.Value,2), round(event.Value,2));
% 
% % Snap sliders to 0.01 THz after release
% fSlider.ValueChangedFcn = @(src,~) ...
%     snapAndRefresh(src, 0.01, fSlider, bwSlider, refresh);
% 
% bwSlider.ValueChangedFcn = @(src,~) ...
%     snapAndRefresh(src, 0.01, fSlider, bwSlider, refresh);
% 
% % Dropdown callbacks
% normMenu.ValueChangedFcn = @(~,~) ...
%     refresh(fSlider.Value, bwSlider.Value);
% 
% compMenu.ValueChangedFcn = @(~,~) ...
%     refresh(fSlider.Value, bwSlider.Value);
% 
% % Initial plot
% refresh(fSlider.Value, bwSlider.Value);


%% 10. Relaxation time vs frequency

% Settings
normalization = 'No pump';  % 'Reference' or 'No pump'

bandwidths = [0 0.1 0.2 0.5];  % THz
minR2 = 0.8;                   % Minimum fit quality
fitStart = 0;                   % Minimum delay for fitting (ps)

% Select normalization
switch normalization
    case 'Reference'
        T = cat(2, Href{idxPump});
    case 'No pump'
        T = cat(2, HnoPump{idxPump});
end

% Frequency and delay arrays
f = freq(:);
t = delayPump(:);

nf = numel(f);
nb = numel(bandwidths);

% Output matrices: frequency x bandwidth
tauAmp = NaN(nf, nb);
tauPhase = NaN(nf, nb);

R2Amp = NaN(nf, nb);
R2Phase = NaN(nf, nb);

% Fit range
fitMask = t >= fitStart;

% Calculate tau(f)
for b = 1:nb

    bw = bandwidths(b);

    for k = 1:nf

        f0 = f(k);

        % Frequency averaging window
        mask = abs(f - f0) <= bw/2;

        if ~any(mask)
            [~, nearest] = min(abs(f-f0));
            mask(nearest) = true;
        end

        % Amplitude averaged over frequency
        amp = mean(abs(T(mask,:)), 1).';

        % Phase of complex-averaged transmission
        phase = unwrap(angle(mean(T(mask,:), 1))).';

        % Fit amplitude
        [tauAmp(k,b), R2Amp(k,b)] = ...
            fitDecay(t(fitMask), amp(fitMask), minR2);

        % Fit phase
        [tauPhase(k,b), R2Phase(k,b)] = ...
            fitDecay(t(fitMask), phase(fitMask), minR2);
    end
end

%% Plot results

figRelaxation = figure('Name', 'Relaxation time vs frequency', ...
       'Color', 'w');

tiledlayout(2,1, 'TileSpacing', 'compact');

colors = lines(nb);

% Amplitude
nexttile;
hold on;

for b = 1:nb
    valid = isfinite(tauAmp(:,b));
    plot(f(valid), tauAmp(valid,b), 'o-', ...
        'Color', colors(b,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('BW = %g THz', bandwidths(b)));
end

ylim([0,4]);
xlabel('Frequency (THz)');
ylabel('\tau (ps)');
title('Amplitude relaxation time');
legend('Location', 'best');
grid on;

% Phase
nexttile;
hold on;

for b = 1:nb
    valid = isfinite(tauPhase(:,b));
    plot(f(valid), tauPhase(valid,b), 'o-', ...
        'Color', colors(b,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('BW = %g THz', bandwidths(b)));
end

ylim([0,4]);
xlabel('Frequency (THz)');
ylabel('\tau (ps)');
title('Phase relaxation time');
legend('Location', 'best');
grid on;

sgtitle(sprintf('Normalization: %s', normalization));

%% Save two relaxation datasets and one combined image

bwTags = arrayfun(@(value) sprintf('%g', value), bandwidths, ...
    'UniformOutput', false);
normTag = regexprep(normalization, '\s+', '_');
relaxationTag = sprintf('%s_norm_%s_BW_%sTHz_minR2_%g_fitStart_%gps', ...
    overviewTag, normTag, strjoin(bwTags, '-'), minR2, fitStart);

if flag_save_csv
    ampColumns = {f};
    phaseColumns = {f};
    fitHeaders = {'frequency_THz'};

    for b = 1:nb
        ampColumns(end+1:end+2) = {tauAmp(:,b), R2Amp(:,b)};
        phaseColumns(end+1:end+2) = {tauPhase(:,b), R2Phase(:,b)};
        fitHeaders(end+1:end+2) = { ...
            sprintf('tau_ps_BW_%gTHz', bandwidths(b)), ...
            sprintf('R2_BW_%gTHz', bandwidths(b))};
    end

    saveColumnsCsv(fullfile(saveDir, ...
        ['09_relaxation_amplitude_' relaxationTag '.csv']), ampColumns, fitHeaders);
    saveColumnsCsv(fullfile(saveDir, ...
        ['10_relaxation_phase_' relaxationTag '.csv']), phaseColumns, fitHeaders);
end

if flag_save_png
    exportgraphics(figRelaxation, ...
        fullfile(saveDir, ['relaxation_' relaxationTag '.png']), ...
        'Resolution', pngResolution);
end


function saveColumnsCsv(filePath, columns, headers)
    % Write the exact labels (including decimal delays and bandwidths) to CSV.
    % Paired time/signal columns may have different lengths.
    nRows = max(cellfun(@numel, columns));
    values = NaN(nRows, numel(columns));

    for k = 1:numel(columns)
        values(1:numel(columns{k}), k) = columns{k}(:);
    end

    writecell([headers; num2cell(values)], filePath, 'Delimiter', ',');
end


function [tau, R2] = fitDecay(t, y, minR2)

    tau = NaN;
    R2 = NaN;

    % Remove invalid values
    valid = isfinite(t) & isfinite(y);
    t = t(valid);
    y = y(valid);

    if numel(t) < 4
        return;
    end

    % Sort by delay
    [t, order] = sort(t);
    y = y(order);

    % Normalize data for numerical stability
    scale = max(y) - min(y);

    if scale < 1e-10 * max(1, max(abs(y)))
        return;
    end

    y0 = y(end);
    yn = (y - y0) / scale;

    % Initial guesses
    tau0 = (max(t) - min(t))/3;
    A0 = yn(1);

    p0 = [A0, log(tau0), 0];

    % Model: A*exp(-t/tau)+B
    model = @(p,x) ...
        p(1)*exp(-x./exp(p(2))) + p(3);

    objective = @(p) ...
        sum((model(p,t-t(1)) - yn).^2);

    options = optimset( ...
        'Display', 'off', ...
        'MaxIter', 2000, ...
        'MaxFunEvals', 5000);

    % Fit
    p = fminsearch(objective, p0, options);

    tauFit = exp(p(2));

    if ~isfinite(tauFit) || tauFit <= 0
        return;
    end

    % Goodness of fit
    yFit = model(p,t-t(1));

    SSE = sum((yn-yFit).^2);
    SST = sum((yn-mean(yn)).^2);

    R2 = 1 - SSE/SST;

    % Reject poor fits
    if isfinite(R2) && R2 >= minR2
        tau = tauFit;
    end

end


function updateTransmissionDelay( ...
    axSpec, axDelay, freq, delay, ...
    Tref, TnoPump, Tbaseline, ...
    f0, bw, normalization, component, ...
    fLabel, bwLabel, tauLabel)

    % Normalization
    switch normalization
        case 'Reference'
            T = Tref;
            base = Tbaseline;
        case 'No pump'
            T = TnoPump;
            base = ones(size(freq));
    end

    % Frequency window
    mask = abs(freq - f0) <= bw/2;

    if ~any(mask)
        [~, k] = min(abs(freq - f0));
        mask(k) = true;
    end

    % Spectral component
    switch component

        case 'Amplitude'
            spectrum = mean(abs(T), 2);
            response = mean(abs(T(mask,:)), 1);
            baseline = mean(abs(base(mask)));
            ylabelText = '|T|';

        case 'Real'
            spectrum = real(mean(T, 2));
            response = real(mean(T(mask,:), 1));
            baseline = real(mean(base(mask)));
            ylabelText = 'Re(T)';

        case 'Imaginary'
            spectrum = imag(mean(T, 2));
            response = imag(mean(T(mask,:), 1));
            baseline = imag(mean(base(mask)));
            ylabelText = 'Im(T)';

        case 'Phase'
            spectrum = angle(mean(T, 2));
            response = angle(mean(T(mask,:), 1));
            baseline = angle(mean(base(mask)));
            ylabelText = 'Phase (rad)';
    end

    %% Spectrum

    cla(axSpec);
    hold(axSpec, 'on');

    plot(axSpec, freq, spectrum, ...
        'k', 'LineWidth', 1.5);

    xline(axSpec, f0-bw/2, 'r--');
    xline(axSpec, f0+bw/2, 'r--');
    xline(axSpec, f0, 'r-');

    xlabel(axSpec, 'Frequency (THz)');
    ylabel(axSpec, ylabelText);
    title(axSpec, 'Transmission spectrum');
    grid(axSpec, 'on');

    %% Response vs delay

    cla(axDelay);
    hold(axDelay, 'on');

    plot(axDelay, delay, response, ...
        'o-', 'LineWidth', 1.5, ...
        'DisplayName', 'Data');

    yline(axDelay, baseline, 'k--', ...
        'DisplayName', 'No pump');

    %% Exponential fit: A*exp(-t/tau)+B

    % Fit only nonnegative delays
    fitMask = delay >= 0 & isfinite(response);

    x = delay(fitMask);
    y = response(fitMask);

    tauLabel.Text = 'Fit unavailable';

    if numel(x) >= 4 && ...
            max(y)-min(y) > 1e-8*max(abs(y))

        % Initial guesses
        A0 = y(1) - y(end);
        B0 = y(end);
        tau0 = (max(x)-min(x))/3;

        % Normalize signal for numerical stability
        scale = max(abs(y));

        % Fit parameters: [A/scale, log(tau), B/scale]
        p0 = [A0/scale, log(tau0), B0/scale];

        model = @(p,t) ...
            p(1)*exp(-t./exp(p(2))) + p(3);

        objective = @(p) ...
            sum((model(p,x)-y/scale).^2);

        options = optimset( ...
            'Display', 'off', ...
            'MaxIter', 2000, ...
            'MaxFunEvals', 5000);

        p = fminsearch(objective, p0, options);

        A = p(1)*scale;
        tau = exp(p(2));
        B = p(3)*scale;

        if isfinite(tau) && tau > 0

            % Smooth fitted curve
            tFit = linspace(min(x), max(x), 300);
            yFit = A*exp(-tFit/tau) + B;

            plot(axDelay, tFit, yFit, 'r-', ...
                'LineWidth', 2, ...
                'DisplayName', 'Exponential fit');

            % Fit quality
            yPred = A*exp(-x/tau) + B;
            SSE = sum((y-yPred).^2);
            SST = sum((y-mean(y)).^2);
            R2 = 1-SSE/SST;

            tauLabel.Text = sprintf( ...
                'tau = %.3f ps     R^2 = %.4f', ...
                tau, R2);

            text(axDelay, 0.04, 0.95, ...
                sprintf('\\tau = %.3f ps', tau), ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'top', ...
                'FontSize', 13, ...
                'FontWeight', 'bold');
        end
    end

    %% Formatting

    xlabel(axDelay, 'Pump-probe delay (ps)');
    ylabel(axDelay, ylabelText);

    title(axDelay, sprintf( ...
        'f = %.2f THz, BW = %.2f THz', f0, bw));

    grid(axDelay, 'on');
    legend(axDelay, 'Location', 'best');

    fLabel.Text = sprintf('f = %.2f THz', f0);
    bwLabel.Text = sprintf('Bandwidth = %.2f THz', bw);

end


function snapAndRefresh(src, step, fSlider, bwSlider, refresh)

    % Snap to exact increments
    value = round(src.Value/step)*step;

    % Respect slider limits
    value = max(src.Limits(1), min(src.Limits(2), value));

    src.Value = value;

    refresh(fSlider.Value, bwSlider.Value);

end


