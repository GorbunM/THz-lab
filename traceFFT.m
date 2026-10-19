function out = traceFFT(A)
    % A = [время, сигнал], равномерная временная сетка
    N = size(A,1);
    dt = mean(diff(A(:,1)));

    S = fft(A(:,2));
    f = ifftshift((-floor(N/2):ceil(N/2)-1).') / (N*dt);
    out = [f S];
end