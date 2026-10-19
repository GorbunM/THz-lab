% folder with data files
curdir = 'C:\Users\mgorbun\Documents\MATLAB\sample10-2\';   

% file with a list of data files
listfile = fullfile(curdir, 'list.txt');

% read list of filenames as string array
A = readmatrix(listfile, 'Delimiter', ',', 'OutputType', 'string');

% frequency range for FFT
frmin=0.15; frmax=2.6;

% create figure and hold for multiple plots
figure
hold on

j = 0;               % counter for valid files
outS = [];           % storage for spectra amplitudes
outPh = [];          % storage for spectra phases

% loop over list of files

for i = A'
    if i{1}(1) ~= '#'                % skip lines starting with '#'
        disp(i{1}(1:2))             % display first 20 chars of filename
        j = j + 1;

        sigfile = fullfile(curdir, i(1));           % build path to current file
        timeTrace = readmatrix(sigfile);    % read data, skip 8 header lines

        specX = TDS_FFT(timeTrace(:, 1:2), frmin, frmax);    % compute spectrum
        
        reffile = fullfile(curdir, i(2));           % build path to current file
        timeTrace = dlmread(reffile, delimeter, 8, 0);    % read data, skip 8 header lines

        ref = TDS_FFT(timeTrace(:, 1:2), frmin, frmax);     % compute spectrum

        freq = spec(:,1);                           % frequency axis
        outS = [outS, spec(:,2)./ref(:,2)];         % accumulate amplitudes
        outPh = [outPh, spec(:,3)];                 % accumulate phases

        % plot(timeTrace(:,1), timeTrace(:,2)-0.02*j, 'LineWidth', 2, 'DisplayName', i{1}(1:end-4));    % plot original time traces with shift
        % plot(freq, (spec(:, 2)./ref(:,2)).^2, 'LineWidth', 1)  % plot spectrum amplitude
        plot(freq, (specX(:, 2)).^2, 'LineWidth', 1, 'DisplayName', i{1}(1:end-4))  % plot power spectrum 
        % plot(freq, spec(:, 3), 'LineWidth', 1)               % plot spectrum phase
    end
end

% ylim([0, 1])
xlabel('THz');
legend('Interpreter','none')

% T = array2table(outS, 'VariableNames', A(:,1)');   % convert accumulated amplitudes to table 
% writetable(T, fullfile(curdir, 'specA.csv'), 'Delimiter', '\t');      % save table as file

% saveas(gcf, fullfile(curdir, 'T.png'))    % save figure to file (disabled)

disp("That's it")