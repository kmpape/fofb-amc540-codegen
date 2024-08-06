function [y_sim, u_sim, x_sim] = sim_lqr_w_constraints(...
            n_samples, n_delay, dist,...
            Ap, Bp, Cp,... % plant
            Ao, Bo, Co, Kf, Kc,... % observer and regulator
            id_to_bpm, id_to_cm,...
            u_max, u_rate,...
            A_awr, B_awr, C_awr, D_awr,...
            ol_mode)
if ~exist('ol_mode','var')
    ol_mode = false;
end

%%
assert((n_delay == 8)||n_delay==9);
use_single = true;
hil_mode = true;

[nx_plant, nu_plant] = size(Bp);
[ny_plant, ~] = size(Cp);
[nx_obs, nu_obs] = size(Bo);
[ny_obs, ~] = size(Co);

% Variables Plant
x_sim_new=zeros(nx_plant,1);
x_sim_old=zeros(nx_plant,1);
y_sim = zeros(ny_plant, n_samples);
u_sim = zeros(nu_plant, n_samples);
x_sim = zeros(nu_plant, n_samples);

% Variables Observer
if use_single == true
    x_pre = single(zeros(nx_obs,1));
    Ao = single(Ao);      Bo = single(Bo);  Co = single(Co);
    Kf = single(Kf);      Kc = single(Kc);
else
    x_pre = double(zeros(nx_obs,1));
    Ao = double(Ao);
    Bo = double(Bo);
    Co = double(Co);
    Kf = double(Kf);
    Kc = double(Kc);
end

% Variables AWR
[ny_awr,nx_awr]=size(C_awr);
x_awr_new=zeros(nx_awr,1);
y_awr=zeros(ny_awr,1);

for k = 1:1:n_samples

    % Measurement
    if ~ol_mode
        y_sim(:, k) = Cp*x_sim_new + dist(:, k);
    else
        y_sim(:, k) = dist(:, k);
    end

    if k > n_delay
        if hil_mode == true
            if use_single == true
                y_meas = single(round(y_sim(id_to_bpm, k-n_delay)*1000,0))*1e-3;
            else
                y_meas = double(round(y_sim(id_to_bpm, k-n_delay)*1000,0))*1e-3;
            end
        else
            if use_single == true
                y_meas = single(y_sim(id_to_bpm, k-n_delay));
            else
                y_meas = double(y_sim(id_to_bpm, k-n_delay));
            end
        end
        
        x_post = x_pre + Kf * (y_meas - Co*x_pre);
        u_obs = -Kc * x_post;
        x_pre = Ao * x_post + Bo * u_obs;
        
        % Compute lower and upper limit
        if use_single == true
            lower_u = max(-u_max,-u_rate+single(y_awr));
            upper_u = min(u_max,u_rate+single(y_awr));
        else
            lower_u = max(-u_max,-u_rate+double(y_awr));
            upper_u = min(u_max,u_rate+double(y_awr));
        end
        
        tmp = max(lower_u, min(upper_u, u_obs));
        u_obs = tmp;
        
        if hil_mode == true
            u_obs = round(u_obs*1e3,0)/1e3;
        end
        u_sim(id_to_cm,k) = double(u_obs);
        
        % AWR
        x_awr_old = x_awr_new;
        x_awr_new = A_awr*x_awr_old + B_awr*(double(tmp));
        y_awr = C_awr*x_awr_new + D_awr*(double(tmp));
    end
    
    % Plant
    x_sim_old = x_sim_new;
    x_sim_new = Ap*x_sim_old + Bp*(u_sim(:,k));
    x_sim(:, k) = x_sim_old;
end
u_sim = u_sim';
y_sim = y_sim';
x_sim = x_sim';

end