function outm = TDS_FFT(timeTrace, frmin, frmax)
    arguments
        timeTrace
        frmin = 0.15
        frmax = 2.5
    end
    
    TT = timeTrace(:, 1:2);
    N = length(TT);
    F = 1 / (TT(2, 1) - TT(1, 1));
    df = F / N;
    
    K  = floor(N/2)+1; 
    freq  = (0:K-1)*df;          % THz

    spec = fft(TT(:, 2));
    if frmax == -1
        frmax = (F-df)/2;
        frmin = df;
    end
    spec = spec(round(frmin/df)+1 : round(frmax/df)+1);
    freq = freq(round(frmin/df)+1 : round(frmax/df)+1);
    outm = [freq', abs(spec), unwrap(angle(spec))]; 
end