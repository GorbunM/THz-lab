curdir = 'C:\Users\mgorbun\Documents\MeasurementKitHost-repair14\Experiments\SpeedTest';

filescsv = dir(fullfile(curdir, '*.csv'));
filesjson = dir(fullfile(curdir, '*.json'));

figure
hold on

% for idx = 1:length(filesjson)
%     j = jsondecode(fileread(fullfile(curdir, filesjson(idx).name)));
%     if j.X == 42.5 && j.Y == 36.9 && j.ScanTimeSeconds == 150 && contains(filesjson(idx).name, 'r0')
%         c = readmatrix(fullfile(curdir, filescsv(idx).name));
%         plot(c(:,1), c(:,2), 'DisplayName', num2str(j.ScanTimeSeconds))
%     end
% end

S = zeros(8000, 1);

for idx = 1:length(filesjson)
    j = jsondecode(fileread(fullfile(curdir, filesjson(idx).name)));
    if j.X == 42.5 && j.Y == 36.9 && j.ScanTimeSeconds == 3 && contains(filesjson(idx).name, 'r0')
        c = readmatrix(fullfile(curdir, filescsv(idx).name));
        filescsv(idx).name
        s = TDS_FFT(c, 0.15, -1);
        S = [S, s(:,2)];
    end
end

s0 = S/length(filesjson);

for idx = 1:length(filesjson)
    j = jsondecode(fileread(fullfile(curdir, filesjson(idx).name)));
    if j.X == 42.5 && j.Y == 36.9 && j.ScanTimeSeconds == 3 && contains(filesjson(idx).name, 'r0')
        c = readmatrix(fullfile(curdir, filescsv(idx).name));
        filescsv(idx).name
        s = TDS_FFT(c, 0.15, -1);
        plot(s(:,1), s(:,2)./s0(:,2), 'DisplayName', num2str(j.ScanTimeSeconds))
    end
end


legend