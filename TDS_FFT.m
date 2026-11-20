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
    % Keep only the non-negative-frequency half of the FFT.  Its last
    % element is the Nyquist frequency, so requested bounds above it must
    % be clipped instead of indexing past the spectrum.
    if frmax == -1
        frmax = (F-df)/2;
        frmin = df;
    end
    firstIndex = max(1, round(frmin/df) + 1);
    lastIndex = min(numel(freq), round(frmax/df) + 1);

    if firstIndex > lastIndex
        error('TDS_FFT:InvalidFrequencyRange', ...
            'Requested range [%.6g, %.6g] THz is outside the available spectrum [%.6g, %.6g] THz.', ...
            frmin, frmax, freq(1), freq(end));
    end

    spec = spec(firstIndex:lastIndex);
    freq = freq(firstIndex:lastIndex);
    outm = [freq', abs(spec), unwrap(angle(spec))]; 
end
