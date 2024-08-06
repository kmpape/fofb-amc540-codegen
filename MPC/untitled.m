nsamples=13;
D1 = randn(1,16);
D2 = randn(1,11);
D3 = randn(1,12);
D4 = randn(1,12);
ybpm = zeros(1,nsamples*3);

nmissing = 0;
kk = 1;
for i = 1:3
    if i == 1
        l_i = length(D1);
        D = D1;
    elseif i==2
        l_i = length(D2);
        D = D2;
    elseif i==3
        l_i = length(D3);
        D = D3;
    end
    
    ybpm(kk:kk+l_i-1) = D;
    if l_i < nsamples
        nmissing = nmissing + nsamples-l_i;
    end
    kk = kk + l_i;
end
if nmissing > 0
    asdf
    ybpm(kk:end) = D4(1,1:nmissing);           
end