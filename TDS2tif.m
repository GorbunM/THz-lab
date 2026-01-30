curdir = 'C:\Users\mgorbun\Documents\MATLAB\TDS\seed\test';
files = dir(fullfile(curdir, 'SLOW X_*.csv'));

X = []; Y = []; vals = struct();
fl = 0
for f = files'
    fl = fl + 1
    name = f.name;
    xy = sscanf(name, 'SLOW X_%f_Y_%f', [1 2]);
    A = dlmread(fullfile(curdir, name), ',', 8, 0);
    spec = TDS_FFT(A(:,1:2));
    spec = spec ./ max(spec);
    
    freq = spec(:,1);
    key = sprintf('X%g_Y%g', xy(1), xy(2));
    key = strrep(key, '.', '_');              
    vals.(key) = spec(:,2);

    X(end+1) = xy(1);
    Y(end+1) = xy(2);
end

xvals = unique(X); yvals = unique(Y);
amp3D = nan(numel(xvals), numel(yvals), numel(freq));

for i = 1:numel(xvals)
    for j = 1:numel(yvals)
        key = sprintf('X%g_Y%g', xvals(i), yvals(j));
        key = strrep(key, '.', '_');
        if isfield(vals, key)
            amp3D(i,j,:) = vals.(key);
        end
    end
end

imwrite(amp3D(:,:,1), 'cube.tif')
for k = 2:size(amp3D,3)
    imwrite(amp3D(:,:,k), 'cube.tif', 'WriteMode', 'append')
end

disp('Thats it')