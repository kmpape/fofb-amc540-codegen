function cols = make_cols(n)
    cols = {'x'};
    for i=1:n
        cols{end+1}=num2str(i);
    end
end