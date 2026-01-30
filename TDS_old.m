% folder with data files
curdir = 'C:\Users\mgorbun\Documents\MATLAB\TDS\seed';   
% file with a list of data files
listfile = fullfile(curdir, 'list.csv');

A = readmatrix(listfile, 'Delimiter', ',', 'OutputType', 'string');

%frequancy sets
frmin=0.15; frmax=1.6; nptfr=146;
freq = frmin:(frmax-frmin)/(nptfr-1):frmax;

figure
hold on

j = 0;

outMatrix = freq';
headers = ['freq'];
for i = A'
    disp(i{1}(1))
    if i{1}(1) ~= '#'
        % disp(i{1}(1))
        j = j + 1;

        sigfile = fullfile(curdir, i(1));
        % timeTrace = dlmread(sigfile, ',', 8, 0);
        % plot(timeTrace(:,1), timeTrace(:,2))
        
        % reffile = fullfile(curdir, i(2));
        % spectrum = trapezFourierT(sigfile, reffile, freq);          % TRANSMITTANCE (choose one ) 
        % spectrum = trapezFourierN(sigfile, reffile, freq, i(4));  % REFRACTIVE INDEX (choose one) 
        spectrum = trapezFourierP(sigfile, freq);              % POWER (choose one)
        outMatrix = [outMatrix, spectrum];
        % headers = [headers, i(1)];
        % plot(freq, spectrum)
    end
end

xlabel('THz');

% T = array2table(outMatrix, 'VariableNames', headers);
% writetable(T, 'C8.csv', 'Delimiter', '\t'); % name the output file

% saveas(gcf, fullfile(curdir, 'T.png'))

% legend

disp("That's it")

