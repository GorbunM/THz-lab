curdir = 'C:\Users\mgorbun\Documents\MATLAB\TDS\seed\test';
files = dir(fullfile(curdir, 'SLOW X_*.csv'));

X = []; 
Y = []; 
Mvals = [];

for f = files'
    name = f.name;
    xy = sscanf(name, 'SLOW X_%f_Y_%f', [1 2]);
    A = dlmread(fullfile(curdir, name), ',', 8, 0);

    X(end+1) = xy(1);
    Y(end+1) = xy(2);
    Mvals(end+1) = max(A(:,2));
end

xvals = unique(X);
yvals = unique(Y);

amp2D = zeros(numel(xvals), numel(yvals));

for n = 1:numel(Mvals)
    ix = find(xvals == X(n));
    iy = find(yvals == Y(n));
    amp2D(ix, iy) = Mvals(n);
end

% ===== PLOT =====
figure
imagesc(xvals, yvals, amp2D.')
axis xy
axis equal tight
colorbar
xlabel('X')
ylabel('Y')
title('Max amplitude map')
colormap(jet)