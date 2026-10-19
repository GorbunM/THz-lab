nobm = readmatrix('C:\\Users\mgorbun\Documents\MATLAB\beamsplitter\a.csv\GVg_-5t5__colimated-no-bs1.dat');
trbm = readmatrix('C:\\Users\mgorbun\Documents\MATLAB\beamsplitter\a.csv\GVg_-5t5__colimated-tr-bs0.dat');
reflbm = readmatrix('C:\\Users\mgorbun\Documents\MATLAB\beamsplitter\a.csv\GVg_-5t5__colimated-refl-bs0.dat');
beforebm = readmatrix('C:\\Users\mgorbun\Documents\MATLAB\beamsplitter\a.csv\GVg_-5t5__colimated-before-bs1.dat');

nobm = nobm(2758:5758, :);
trbm = trbm(203:3203, :);
reflbm = reflbm(1042:4042, :);
beforebm = beforebm(645:3645, :);

iter = {nobm, trbm, reflbm, beforebm};
names = {'nobm', 'trbm', 'reflbm', 'beforebm'};

for idx = 1:numel(iter)
    T = iter{idx};
    if T(1,7) == 310
        iter{idx} = flip(T, 1);
    end
end
%%
figure
hold on

for idx = 1:numel(iter)
    T = iter{idx};
    plot(T(:, 7), T(:, 1), 'DisplayName', names{idx})
end
legend
%%
% figure
% hold on
% 
% for idx = 1:numel(iter)
%     T = iter{idx};
%     plot(T(:, 7), T(:, 1)./beforebm(:, 1), 'DisplayName', names{idx})
% end
% legend
% ylim([0, 1])
%%
figure
hold on
for idx = 1:2
    T = iter{idx+1};
    plot(T(:, 7), T(:, 1)./nobm(:, 1), 'DisplayName', names{idx+1})
    mean(T(:, 1)./nobm(:, 1))
end
legend
ylim([0, 1])
%%

figure
hold on
for idx = 2
    T = iter{idx};
    plot(T(:, 7), T(:, 1)./reflbm(:, 1), 'DisplayName', names{idx})
    mean(T(:, 1)./nobm(:, 1))
end
legend
ylim([0, 1])