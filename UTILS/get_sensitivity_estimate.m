function [fs, S] = get_sensitivity_estimate(Y, D, nfft, Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the Welch method to obtain an estimate of the sensitivity Y./D      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [N, ny] = size(Y);
    assert(N>=nfft);
    navg = floor(N/nfft);
    numfreqs = floor(nfft / 2 + 1);
    S = zeros(navg, numfreqs );
    for i=1:navg
        Yi = Y(1+(i-1)*nfft:i*nfft,:);
        Di = D(1+(i-1)*nfft:i*nfft,:);
        
        [fs, Ypsd, ~] = get_all_psd(Yi, nfft, Fs);
        [~,  Dpsd, ~] = get_all_psd(Di, nfft, Fs);
        
        Ypsdnorm = vecnorm(Ypsd,2,1);
        Dpsdnorm = vecnorm(Dpsd,2,1);
        
        S(i,:) = sqrt(Ypsdnorm ./ Dpsdnorm);
        
    end

    S = mean(S,1);
end



