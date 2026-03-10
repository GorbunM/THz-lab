c = 2.99792458e8;

fmin = 1e3;
fmax = 1e22;

figure

%% верхняя шкала — частота
ax1 = axes;
set(ax1,'XScale','log','XAxisLocation','top')
xlim([fmin fmax])
yticks([])
xlabel('Frequency (Hz)')
grid on
pos = ax1.Position;

%% нижняя шкала — длина волны
ax2 = axes;
set(ax2,'Position',pos,...
        'Color','none',...
        'XScale','log',...
        'XAxisLocation','bottom',...
        'YTick',[],...
        'XDir','reverse')
xlim([fmin fmax])
xlabel('Wavelength (m)')

%% ось для диапазонов
ax3 = axes;
set(ax3,'Position',pos,...
        'Color','none',...
        'XScale','log',...
        'XTick',[],...
        'YTick',[],...
        'Box','off')
xlim([fmin fmax])
ylim([0 1])
hold on

bands = {
'Radio',      1e3, 3e9
'Microwave',  3e9, 3e11
'THz',        3e11, 3e13
'IR',         3e13, 4e14
'Visible',    4e14, 7.5e14
'UV',         7.5e14, 3e16
'X-ray',      3e16, 3e19
'Gamma',      3e19, 1e22
};

for i = 1:size(bands,1)

    x1 = bands{i,2};
    x2 = bands{i,3};

    patch([x1 x2 x2 x1],[0 0 1 1],...
          [0.85 0.85 0.85],'EdgeColor','k')

    text(sqrt(x1*x2),0.5,bands{i,1},...
        'HorizontalAlignment','center')
end

%% тики частоты (10^n)
fticks = 10.^(3:3:21);
set(ax1,'XTick',fticks)

%% тики длины волны (10^n m)
lambda_ticks = 10.^(-12:3:6);
fticks_lambda = c ./ lambda_ticks;

set(ax2,'XTick',fticks_lambda)
set(ax2,'XTickLabel',compose('10^{%d}',-12:3:6))