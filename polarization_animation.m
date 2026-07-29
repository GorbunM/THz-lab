phi = linspace(0, 10*pi, 1000);
phase = 0.25*pi;
mX = 1;
mY = 2;

%%

v = VideoWriter('R.mp4','MPEG-4');
v.FrameRate = 20;
open(v)

figure
for k = 1:numel(phi)
    hold on
    Phi = phi(1:k);
    z = mX*cos(Phi) + 1i*mY*sin(Phi+phase);
    magn = max(mX, mY);
    plot(z, 'r'); 
    z0 = z(k);
    compass(z0)
    axis equal
    ylim([-magn, magn]); 
    xlim([-magn, magn]);
    frame = getframe(gcf);
    writeVideo(v, frame);
    drawnow
    hold off
    clf
end

close(v)
disp('That''s it')
