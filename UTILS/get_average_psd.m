function [fs, psd, cs] = get_average_psd(Y, nfft, F_S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a matrix Y dim(n_samples, n_measurements), it returns a psd with  %
% the maximum amplitute over all measurements.                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, n_y] = size(Y);

    numfreqs = floor(nfft / 2 + 1);
    Psd = zeros(n_y, numfreqs);
    Cs = zeros(n_y, numfreqs -1);
    % Psd = zeros(n_y, numfreqs);

    for iy = 1 : n_y
        [fs, Psd(iy, :), Cs(iy, :)] = mypsd(Y(:,iy), nfft, F_S);
    end

    psd = mean(Psd, 1);
    cs = mean(Cs, 1);
end



