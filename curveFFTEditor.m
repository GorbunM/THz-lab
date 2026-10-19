function [x, y, f, F] = curveFFTEditor()
    % curveFFTEditor()
    %
    % Верх: интерактивное редактирование кривой
    % Низ: правая половина abs(FFT)
    %
    % Управление:
    %   мышь  — тянуть кривую вверх/вниз
    %   Enter — закончить и вернуть x, y, f, F
    %
    % Выход:
    %   x — координата
    %   y — итоговая кривая
    %   f — частотная ось5
    %   F — комплексный односторонний FFT

    N = 512;

    x = linspace(0, 1, N);
    dx = x(2) - x(1);

    y = zeros(1, N);

    sigma = 0.03;      % ширина утягивания соседей
    yStartLimit = 1;   % только начальный масштаб

    fig = figure( ...
        'Name', 'Curve + FFT editor', ...
        'NumberTitle', 'off');

    tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');

    % ============================================================
    % Верхний график: редактируемая кривая
    % ============================================================

    ax1 = nexttile;
    hold(ax1, 'on');
    grid(ax1, 'on');
    box(ax1, 'on');

    hLine = plot(ax1, x, y, 'LineWidth', 2);
    hPts  = plot(ax1, x, y, '.', 'MarkerSize', 8);

    xlim(ax1, [x(1), x(end)]);
    ylim(ax1, [-yStartLimit, yStartLimit]);

    ylabel(ax1, 'y(x)');
    title(ax1, 'Drag curve vertically. Press Enter when done.');

    % ============================================================
    % Нижний график: правая половина спектра
    % ============================================================

    ax2 = nexttile;
    hold(ax2, 'on');
    grid(ax2, 'on');
    box(ax2, 'on');

    [f, F, A] = calcFFT(y, dx);

    hFFT = plot(ax2, f, A, 'LineWidth', 2);

    xlim(ax2, [0, max(f)]);
    ylim(ax2, [0, 1]);

    xlabel(ax2, 'frequency');
    ylabel(ax2, '|FFT|');
    title(ax2, 'Right half of Fourier spectrum');

    % ============================================================
    % Состояние мыши
    % ============================================================

    dragging = false;
    lastY = 0;
    idx0 = 1;

    fig.WindowButtonDownFcn   = @mouseDown;
    fig.WindowButtonUpFcn     = @mouseUp;
    fig.WindowButtonMotionFcn = @mouseMove;
    fig.KeyPressFcn           = @keyPress;

    updateFFT();

    uiwait(fig);

    if isvalid(fig)
        close(fig);
    end

    % ============================================================
    % Callback-функции
    % ============================================================

    function mouseDown(~, ~)
        cp = ax1.CurrentPoint;

        mx = cp(1, 1);
        my = cp(1, 2);

        if mx < x(1) || mx > x(end)
            return
        end

        [~, idx0] = min(abs(x - mx));

        dragging = true;
        lastY = my;
    end

    function mouseMove(~, ~)
        if ~dragging
            return
        end

        cp = ax1.CurrentPoint;
        my = cp(1, 2);

        dy = my - lastY;
        lastY = my;

        % Вес влияния на соседние точки
        w = exp(-((x - x(idx0)).^2) / (2*sigma^2));

        % Без ограничения по амплитуде
        y = y + dy*w;

        hLine.YData = y;
        hPts.YData = y;

        expandYLim(ax1, y);

        updateFFT();

        drawnow limitrate
    end

    function mouseUp(~, ~)
        dragging = false;
    end

    function keyPress(~, event)
        if strcmp(event.Key, 'return')
            uiresume(fig);
        end
    end

    % ============================================================
    % Обновление FFT
    % ============================================================

    function updateFFT()
        [f, F, A] = calcFFT(y, dx);

        hFFT.XData = f;
        hFFT.YData = A;

        ymax = max(A);

        if ymax <= 0
            ymax = 1;
        end

        ylim(ax2, [0, 1.05*ymax]);
    end

    % ============================================================
    % Авторасширение вертикальной оси
    % ============================================================

    function expandYLim(ax, ydata)
        ymin = min(ydata);
        ymax = max(ydata);

        lim = ylim(ax);

        if ymin < lim(1) || ymax > lim(2)
            pad = 0.1 * max(ymax - ymin, 1);

            newLim = [
                min(lim(1), ymin - pad), ...
                max(lim(2), ymax + pad)
            ];

            ylim(ax, newLim);
        end
    end
end

% ================================================================
% Односторонний FFT
% ================================================================

function [f, F, A] = calcFFT(y, dx)
    N = numel(y);

    Ffull = fft(y);

    Nhalf = floor(N/4) + 1;

    F = Ffull(1:Nhalf);

    A = abs(F) / N;

    % Амплитудная поправка для одностороннего спектра
    if Nhalf > 2
        A(2:end-1) = 2*A(2:end-1);
    end

    f = (0:Nhalf-1) / (N*dx);
end