f = @(x) parametricTMM(x);   
X = linspace(0.1,2.5,100); 
%%

v = VideoWriter('A.mp4','MPEG-4');
v.FrameRate = 5;
open(v)

figure
for k = 1:numel(X)
    [nu, R, eps, A] = f(X(k));
    [ymax, idx] = findpeaks(A);
    yfit = interp1(nu(idx), ymax, nu, 'spline');
    tiledlayout(2,1)
    nexttile, plot(nu, real(eps), 'r', nu, imag(eps), 'b'); ylim([-40, 75])
    title(sprintf('omega_p = 2pi x %.2f', floor(X(k)*10)))
    nexttile, plot(nu, A, '-', nu, yfit, '-'); ylim([0, 1])
    frame = getframe(gcf);
    writeVideo(v, frame);
    drawnow
end

close(v)
disp('That''s it')