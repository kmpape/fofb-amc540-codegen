function [Sfiltered, Sorig, w_Hz] = estim_sensitivity(yon, doff, n_samples, Fs)
    assert(size(yon, 2) == n_samples);
    assert(size(doff, 2) == n_samples);
    yon = yon - mean(yon, 2);
    doff = doff - mean(doff, 2);
    [Sorig,F] = tfestimate(doff',yon', n_samples, [], n_samples);
    w_Hz = F*Fs/2/pi;
    ind_csv = w_Hz <= Fs/2; 
    w_Hz = w_Hz(ind_csv);
    n_csv = sum(ind_csv);
    Sfiltered = zeros(n_csv, size(Sorig,2));
    for ii=1:size(Sfiltered,2)
        tmp = sgolayfilt(20*log10(abs(Sorig(:,ii))), 5, 501);
        Sfiltered(:,ii) = tmp(ind_csv);
    end
end
