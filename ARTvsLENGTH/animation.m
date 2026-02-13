nunu = nu*1e-12;
x_layers = cumsum(d_layers) * 1e4;

col = {'y', 'c', 'b'}
names = {'Au', 'Si', 'XXX'}
%%

v = VideoWriter('R.mp4','MPEG-4');
v.FrameRate = 10;
open(v)

figure
for k = 1:numel(thcEXP)
    tiledlayout(2,1)
    nexttile, plot(nunu, refl(:,k), 'r'); ylim([0, 1])
    % title(sprintf('omega_p = 2pi x %.2f', floor(X(k)*10)))    
    title(sprintf('d = %.2f', thcEXP(k)*1e4))
    nexttile
    x = x_layers(:, k);
    for i = 1:3
        p = patch([x(i) x(i+1) x(i+1) x(i)],[0 0 1 1],col{i}, 'DisplayName', names{i}); ylim([0, 1]); xlim([0, x_layers(end)])
        % hatchfill2(p, diagonal)
    end
    legend
    frame = getframe(gcf);
    writeVideo(v, frame);
    drawnow
end

close(v)
disp('That''s it')

%%
% 
% figure
% tiledlayout(2,1)
% nexttile, plot(nu, refl(:,1), 'r'); ylim([0, 1])
% nexttile, plot(nu, trans(:,1), 'b'); ylim([0, 1])

%%

% figure
% plot(real(eps_layers(2,:)), 'k--')
% % plot(imag(eps_layers(2,:)), 'b-')