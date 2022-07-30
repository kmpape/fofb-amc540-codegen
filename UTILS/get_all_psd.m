function [fs, psd, cs] = get_all_psd(Y, nfft, F_S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a matrix Y dim(n_samples, n_measurements), it returns a psd with  %
% the maximum amplitute over all measurements.                            %
% TIME AXIS = 1 (rows)                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, n_y] = size(Y);

    numfreqs = floor(nfft / 2 + 1);
    psd = zeros(n_y, numfreqs );
    cs = zeros(n_y, numfreqs-1 );
    % Psd = zeros(n_y, numfreqs);

    for iy = 1 : n_y
        [fs, psd(iy, :), cs(iy, :)] = mypsd(Y(:,iy), nfft, F_S);
    end
end



