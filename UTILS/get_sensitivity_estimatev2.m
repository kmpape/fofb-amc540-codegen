function [fs, S] = get_sensitivity_estimatev2(Y, D, nfft, Fs, xpad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the Welch method to obtain an estimate of the sensitivity Y./D      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [N, ny] = size(Y);
    assert(N>=nfft);
    navg = floor(N/nfft);
    numfreqs = floor((nfft*xpad) / 2 + 1);    
    delta_f = Fs / (nfft*xpad);
    fs = delta_f * (1:numfreqs-1);
    S = zeros(navg, numfreqs-1);
    for i=1:navg
        Yi = Y(1+(i-1)*nfft:i*nfft,:).*hamming(nfft);
        Di = D(1+(i-1)*nfft:i*nfft,:).*hamming(nfft);
        
        Ypsd = 2/(nfft*xpad)*abs(fft(Yi, xpad*nfft, 1)).^2;
        Ypsd = Ypsd(2:numfreqs,:)';
        
        Dpsd = 2/(nfft*xpad)*abs(fft(Di, xpad*nfft, 1)).^2;
        Dpsd = Dpsd(2:numfreqs,:)';
        
        Ypsdnorm = vecnorm(Ypsd,2,1);
        Dpsdnorm = vecnorm(Dpsd,2,1);
        
        S(i,:) = sqrt(Ypsdnorm ./ Dpsdnorm);
        
    end

    S = mean(S,1);
end



