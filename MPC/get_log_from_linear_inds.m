function inds = get_log_from_linear_inds(w_log,w_lin)
    ilin = 0;
    n_log = length(w_log);
    n_lin = length(w_lin);
    inds = zeros(n_log,1);
    for ilog=1:n_log
        while (1)
            ilin = ilin + 1;
            if w_lin(ilin) >= w_log(ilog)
                inds(ilog) = ilin;
                break;
            end            
        end
    end
end