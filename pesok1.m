f = @(x) parametricTMM(x);   
X = [0.37, 0.97, 2.37]; 
%%


for k = 1:numel(X)
    figure
    [nu, R, eps, A] = f(X(k));
    [ymax, idx] = findpeaks(A);
    yfit = interp1(nu(idx), ymax, nu, 'spline');
    tiledlayout(2,1)
    nexttile, plot(nu, real(eps), 'r', nu, imag(eps), 'b'); ylim([-40, 75])
    % title(sprintf('omega_p = 2pi x %.2f', floor(X(k)*10)))    
    title(sprintf('x = %.2f', X(k)))
    nexttile, plot(nu, A, '-', nu, yfit, '-'); ylim([0, 1])
end

close(v)
disp('That''s it')