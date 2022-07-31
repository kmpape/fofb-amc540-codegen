function [fs, psd, cs] = mypsd(y, nfft, F_S)
% function [fs, psd, cs] = mypsd(y, nfft, F_S)
% 
% THIS IS THE ONLY WAY TO DO IT!
% tested against all known plots and tools
%
% fs  : frequencies
% psd : psd in [units^2/Hz]
% cs  : integrated spectrum [units]
% 
% y   : input [units] 
% nfft: FFT block size
% F_S : sample frequency [Hz]
%
if nargin < 4
    take_cumsum = true;
end

numfreqs = floor(nfft / 2 + 1);
N = length(y);
% step could be different if we use overlap and window
step = nfft;
start = 1;
n = 1;

psd = zeros(floor(N / nfft), nfft);

while 1
    last = start + nfft - 1;
    if last > N
        break
    end
    block = y(start:last);
    pp = 2 / length(block) * abs(fft(block)).^2;
    psd(n, :) = pp;
    start = start + step;
    n = n + 1;
end

% average blocks
psd = mean(psd, 1);

delta_f = F_S / nfft;

fs = delta_f * (0:numfreqs-1);

% power per Hz
psd = psd(1:numfreqs) / F_S;

% remove 0 frequency % CHANGED for psd
fs = fs(1:end);
psd = psd(1:end);
% integrate: limit as df -> 0 of sum(psd * delta_f)
cs = sqrt(cumsum(psd(2:end) * delta_f));

