f = @(x) parametricTMM(x);   
X = linspace(1, 2000, 200); 
% X = 1;
%%
delete('R.mp4')
v = VideoWriter('R.mp4','MPEG-4');
v.FrameRate = 35;
open(v)

figure
for k = 1:numel(X)
    [nu, R, eps, A, T, phi] = f(X(k));
    % [ymax, idx] = findpeaks(A);
    % yfit = interp1(nu(idx), ymax, nu, 'spline');
    % tiledlayout(2,1)
    % nexttile, plot(nu, real(eps), 'r', nu, imag(eps), 'b'); ylim([-100, 100])
    % title(sprintf('omega_p = 2pi x %.2f', floor(X(k)*10)))    
    % title(sprintf('x = %.2f', X(k)))
    % title(sprintf('n_{Di} = %.2f', 3.5/X(k)))
    % nexttile, plot(nu, R, '-', nu, yfit, '-'); ylim([0, 1])
    % nexttile, hold on;
    clf
    hold on;
        plot(nu, R, 'r-', 'DisplayName', 'R', 'LineWidth', 2); 
        plot(nu, T, 'b-', 'DisplayName', 'T', 'LineWidth', 2); 
        plot(nu, A, 'k-', 'DisplayName', 'A', 'LineWidth', 2);
        ylim([0, 1]); legend; hold off
    ax = gca;
    ax.FontSize = 24;
    % nexttile, plot(nu, A, 'k-', 'DisplayName', 'A'); ylim([0, 1]); legend
    % plot(nu, phi, nu, T); ylim([-pi, pi])
    frame = getframe(gcf);
    writeVideo(v, frame);
    drawnow
end

close(v)

% saveas(gcf, '3.png')

disp('That''s it')