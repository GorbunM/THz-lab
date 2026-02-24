f = @(x) parametricTMM(x);   
X = linspace(1, 7/3, 100); 
%%

v = VideoWriter('R.mp4','MPEG-4');
v.FrameRate = 20;
open(v)

figure
for k = 1:numel(X)
    [nu, R, eps, A, T] = f(X(k));
    % [ymax, idx] = findpeaks(A);
    % yfit = interp1(nu(idx), ymax, nu, 'spline');
    tiledlayout(2,1)
    nexttile, plot(nu, real(eps), 'r', nu, imag(eps), 'b'); ylim([-10, 10])
    % title(sprintf('omega_p = 2pi x %.2f', floor(X(k)*10)))    
    % title(sprintf('x = %.2f', X(k)))
    title(sprintf('n_{Di} = %.2f', 3.5/X(k)))
    % nexttile, plot(nu, R, '-', nu, yfit, '-'); ylim([0, 1])
    nexttile, plot(nu, R, '-', nu, T, '-'); ylim([0, 1])
    frame = getframe(gcf);
    writeVideo(v, frame);
    drawnow
end

close(v)
disp('That''s it')