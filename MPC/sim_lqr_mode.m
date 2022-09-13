function [y_sim, u_sim] = sim_lqr_mode(...
            n_samples, n_delay, dist,...
            Ap, Bp, Cp,... % plant
            Aoi, Boi, Coi, Kfi, Kci, UR, VR,... % observer and regulator
            id_to_bpm, id_to_cm,...
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
[nx_obs, nu_obs, ~] = size(Boi);
[~, ~, ny_obs] = size(Coi);

% Variables Plant
x_sim_new=zeros(nx_plant,1); 
x_sim_old=zeros(nx_plant,1);
y_sim = zeros(ny_plant, n_samples);
u_sim = zeros(nu_plant, n_samples);

% Variables Observer
if use_single == true
    x_pre = single(zeros(nx_obs,1,ny_obs));
    u_mode = single(zeros(nu_obs,1));
    Aoi = single(Aoi);	Boi = single(Boi);  Coi = single(Coi);
    Kfi = single(Kfi);	Kci = single(Kci);
else
    x_pre = double(zeros(nx_obs,1,ny_obs));
    u_mode = double(zeros(nu_obs,1));
    Aoi = double(Aoi);	Boi = double(Boi);  Coi = double(Coi);
    Kfi = double(Kfi);  Kci = double(Kci);
end

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
        y_mode = UR'*y_meas;
        for imode = 1:ny_obs
            x_post = x_pre(:,:,imode) + Kfi(:,:,imode) * (y_mode(imode) - Coi(:,:,imode)*x_pre(:,:,imode));
            u_mode(imode,1) = -Kci(:,:,imode) * x_post;
            x_pre(:,:,imode) = Aoi(:,:,imode) * x_post + Boi(:,:,imode) * u_mode(imode);
        end
        u_obs = VR*u_mode;

        if hil_mode == true
            u_obs = round(u_obs*1e6,0)/1e6;
        end
        u_sim(id_to_cm,k) = double(u_obs);
    end
    
    % Plant
    x_sim_old = x_sim_new;
    x_sim_new = Ap*x_sim_old + Bp*(u_sim(:,k));
end
u_sim = u_sim';
y_sim = y_sim';

end