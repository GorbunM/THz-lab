curdir = 'C:\Users\mgorbun\Documents\Presentations\';
figDir = 'fig';
outFile = fullfile(curdir, 'slides1.md');

ext = {'*.png','*.jpg','*.jpeg','*.tif','*.tiff','*.bmp'};

files = [];
for k = 1:numel(ext)
    files = [files; dir(fullfile(curdir, figDir, ext{k}))];
end

fid = fopen(outFile,'w');

for k = 1:numel(files)
    name = files(k).name;
    % fprintf(fid,'#\n');
    fprintf(fid,'# %s\n\n', name(1:3));
    % fprintf(fid,'---\n');    
    % fprintf(fid,'layout: blank \n');
    % fprintf(fid,'---\n');
    fprintf(fid,'![](%s/%s){width=6cm}\n\n', figDir, name);

    if k < numel(files)
        fprintf(fid,'---\n\n');
    end
end

fclose(fid);


%%

curdir = 'C:\Users\mgorbun\Documents\Presentations\';
figDir = 'fig';
outFile = fullfile(curdir, 'slides1.org');

ext = {'*.png','*.jpg','*.jpeg'}%,'*.tif','*.tiff','*.bmp'};

files = [];
for k = 1:numel(ext)
    files = [files; dir(fullfile(curdir, figDir, ext{k}))];
end

fid = fopen(outFile,'w');

fprintf(fid,'#+TITLE: Slides\n\n');

for k = 1:numel(files)
    name = files(k).name;

    title = name(1:min(3,numel(name)));

    % fprintf(fid,'* %s\n\n', title);
    fprintf(fid,'#+ATTR_PANDOC: :width 16cm\n');
    fprintf(fid,'[[file:%s/%s]]\n\n', figDir, name);
end

fclose(fid);
disp('That''s it')

%%

curdir = 'C:\Users\mgorbun\Documents\Presentations\fig';

files = dir(fullfile(curdir,'*.tif'));

for k = 1:numel(files)
    name = files(k).name;

    I = imread(fullfile(curdir,name));
    imwrite(I, fullfile(curdir, [name(1:end-4) '.png']));
end

disp('1')