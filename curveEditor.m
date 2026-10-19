function y = curveEditor()
    % curveEditor()
    % Интерактивное вертикальное редактирование кривой.
    % На выходе y — массив значений итоговой кривой.

    N = 200;
    x = linspace(0, 1, N);
    y = zeros(1, N);

    sigma = 0.04;      % ширина "утягивания" соседей
    ampLimit = 1;      % предел по y

    fig = figure( ...
        'Name', 'Curve editor', ...
        'NumberTitle', 'off');

    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    hLine = plot(ax, x, y, 'LineWidth', 2);
    hPts  = plot(ax, x, y, '.', 'MarkerSize', 10);

    ylim(ax, [-ampLimit ampLimit]);
    xlim(ax, [0 1]);

    title(ax, 'Drag curve vertically. Press Enter when done.');

    dragging = false;
    lastY = 0;
    idx0 = 1;

    fig.WindowButtonDownFcn = @mouseDown;
    fig.WindowButtonUpFcn   = @mouseUp;
    fig.WindowButtonMotionFcn = @mouseMove;
    fig.KeyPressFcn = @keyPress;

    uiwait(fig);

    if isvalid(fig)
        close(fig);
    end

    % ---------- callbacks ----------

    function mouseDown(~, ~)
        cp = ax.CurrentPoint;
        mx = cp(1, 1);
        my = cp(1, 2);

        [~, idx0] = min(abs(x - mx));

        dragging = true;
        lastY = my;
    end

    function mouseMove(~, ~)
        if ~dragging
            return
        end

        cp = ax.CurrentPoint;
        my = cp(1, 2);

        dy = my - lastY;
        lastY = my;

        % вес влияния на соседей
        w = exp(-((x - x(idx0)).^2) / (2*sigma^2));

        y = y + dy*w;
        y = max(min(y, ampLimit), -ampLimit);

        hLine.YData = y;
        hPts.YData = y;

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
end