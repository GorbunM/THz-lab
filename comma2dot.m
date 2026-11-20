%% Replace decimal commas with dots

curdir = "C:\Users\mgorbun\Documents\MATLAB\2026-10-30 cnt90 1113 nm repeat with less gap\2026-10-30 cnt90 1113 nm repeat with less gap\perpendicular 1113 nm";

files = dir(fullfile(curdir, '*.dat'));

for i = 1:numel(files)

    filename = fullfile(curdir, files(i).name);

    txt = fileread(filename);

    % Replace commas between digits
    txt = regexprep(txt, '(?<=\d),(?=\d)', '.');

    fid = fopen(filename, 'w');
    fwrite(fid, txt, 'char');
    fclose(fid);

    fprintf('Processed: %s\n', files(i).name);
end