curdir = 'C:\Users\mgorbun\Downloads\Telegram Desktop\auto\auto';
files = dir(fullfile(curdir, "auto_*_*_0*"));
x = []; y = []; z = [];

for k = 1:numel(files)
    ij = sscanf(files(k).name, 'auto_%d_%d_0');
    data = readmatrix(fullfile(curdir, files(k).name));
    x(k) = ij(1);
    y(k) = ij(2);
    z(k) = max(data(:,2));
end

[X,Y] = meshgrid(unique(x), unique(y));
Z = nan(size(X));
for k = 1:numel(z)
    Z(Y == y(k) & X == x(k)) = z(k);
end
imagesc(unique(x), unique(y), Z)
axis xy equal tight
colorbar

%%

