function [y_sim,u_sim,x_sim,...
    obs_y,obs_u,obs_x0,obs_xd,...
    fgm_x0,fgm_xd,fgm_u,fgm_out,lower_u,upper_u] = sim_mpc(...
            n_samples, n_delay, dist,...
            Ap, Bp, Cp,... % Plant
            Ao, Bo, Co, Ad, Cd, LxN_obs, Lxd_obs,... % Observer
            J_MPC, q_mat, beta_fgm,... % FGM
            u_max, u_rate,...
            id_to_bpm, id_to_cm,...
            A_awr, B_awr, C_awr, D_awr,...
            SOFB_setp, ol_mode)
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
x_sim_new=zeros(nx_plant,1); x_sim_old=zeros(nx_plant,1);
y_sim = zeros(ny_plant, n_samples);
u_sim = zeros(nu_plant, n_samples);
x_sim = zeros(nu_plant, n_samples);

% Variables AWR
[ny_awr,nx_awr]=size(C_awr);
x_awr_new=zeros(nx_awr,1);
y_awr=zeros(ny_awr,1);

% Variables Observer
if use_single == true
    x0_obs_old=single(zeros(nx_obs,1)); x1_obs_old=single(zeros(nx_obs,1));
    x2_obs_old=single(zeros(nx_obs,1)); x3_obs_old=single(zeros(nx_obs,1));
    x4_obs_old=single(zeros(nx_obs,1)); x5_obs_old=single(zeros(nx_obs,1));
    x6_obs_old=single(zeros(nx_obs,1)); x7_obs_old=single(zeros(nx_obs,1));
    x8_obs_old=single(zeros(nx_obs,1));
    xd_obs_old=single(zeros(ny_obs,1));
    
    Apow1 = single(Ao.^1);   Apow2 = single(Ao.^2);
    Apow3 = single(Ao.^3);   Apow4 = single(Ao.^4);
    Apow5 = single(Ao.^5);   Apow6 = single(Ao.^6);
    Apow7 = single(Ao.^7);   Apow8 = single(Ao.^8);
    Apow9 = single(Ao.^9);
    Ao = single(Ao);      Bo = single(Bo);  Co = single(Co);
    Ad = single(Ad);      Cd = single(Cd);
    LxN_obs = single(LxN_obs);  Lxd_obs = single(Lxd_obs);
    
    obs_y = single(zeros(ny_obs,n_samples));
    obs_u = single(zeros(nu_obs,n_samples));
    obs_xd = single(zeros(ny_obs,n_samples));
    obs_x0 = single(zeros(nx_obs,n_samples));

    fgm_x0 = single(zeros(nx_obs,n_samples));
    fgm_xd = single(zeros(ny_obs,n_samples));
    fgm_u = single(zeros(nu_obs,n_samples));
    fgm_out = single(zeros(nu_obs,n_samples));
else
    x0_obs_old=double(zeros(nx_obs,1)); x1_obs_old=double(zeros(nx_obs,1));
    x2_obs_old=double(zeros(nx_obs,1)); x3_obs_old=double(zeros(nx_obs,1));
    x4_obs_old=double(zeros(nx_obs,1)); x5_obs_old=double(zeros(nx_obs,1));
    x6_obs_old=double(zeros(nx_obs,1)); x7_obs_old=double(zeros(nx_obs,1));
    x8_obs_old=double(zeros(nx_obs,1));
    xd_obs_old=double(zeros(ny_obs,1));
    
    Apow1 = double(Ao.^1);   Apow2 = double(Ao.^2);
    Apow3 = double(Ao.^3);   Apow4 = double(Ao.^4);
    Apow5 = double(Ao.^5);   Apow6 = double(Ao.^6);
    Apow7 = double(Ao.^7);   Apow8 = double(Ao.^8);
    Apow9 = double(Ao.^9);
    Ao = double(Ao);      Bo = double(Bo);  Co = double(Co);
    Ad = double(Ad);      Cd = double(Cd);
    LxN_obs = double(LxN_obs);  Lxd_obs = double(Lxd_obs);
    
end

% Variables Fast Gradient
if use_single == true
    J = single(J_MPC);
    beta_fgm = single(beta_fgm);
    q_mat = single(q_mat);
    z_new=single(zeros(nu_obs,1));
else
    J = double(J_MPC);
    beta_fgm = double(beta_fgm);
    q_mat = double(q_mat);
    z_new=double(zeros(nu_obs,1));
end

MAX_ITER = 20;
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
        obs_y(:,k-n_delay) = y_meas;
        % Observer - state update
        if use_single == true
            obs_u(:,k-n_delay) = single(u_sim(id_to_cm,k-1));
            x0_obs_new = Ao*x0_obs_old + Bo*single(u_sim(id_to_cm,k-1));
        else
            x0_obs_new = Ao*x0_obs_old + Bo*double(u_sim(id_to_cm,k-1));
        end
        xd_obs_new = Ad*xd_obs_old;
        x1_obs_new = x0_obs_old;
        x2_obs_new = x1_obs_old;
        x3_obs_new = x2_obs_old;
        x4_obs_new = x3_obs_old;
        x5_obs_new = x4_obs_old;
        x6_obs_new = x5_obs_old;
        x7_obs_new = x6_obs_old;
        x8_obs_new = x7_obs_old;
        if (n_delay==9); x9_obs_new = x8_obs_old; end

        % Observer - measurement update
        if (n_delay==9)
            delta_y = y_meas - Co*x9_obs_new - Cd*xd_obs_new;
        else
            delta_y = y_meas - Co*x8_obs_new - Cd*xd_obs_new;
        end
        delta_xN = LxN_obs * delta_y;
        delta_xd = Lxd_obs * delta_y;
        xd_obs_new = xd_obs_new + delta_xd;         % -> result
        if (n_delay==9)
            x9_obs_new = x9_obs_new + delta_xN;      
            x8_obs_new = x8_obs_new + Apow1 * delta_xN;
            x7_obs_new = x7_obs_new + Apow2 * delta_xN;
            x6_obs_new = x6_obs_new + Apow3 * delta_xN;
            x5_obs_new = x5_obs_new + Apow4 * delta_xN;
            x4_obs_new = x4_obs_new + Apow5 * delta_xN;
            x3_obs_new = x3_obs_new + Apow6 * delta_xN;
            x2_obs_new = x2_obs_new + Apow7 * delta_xN;
            x1_obs_new = x1_obs_new + Apow8 * delta_xN; 
            x0_obs_new = x0_obs_new + Apow9 * delta_xN; % -> result
        else
            x8_obs_new = x8_obs_new + delta_xN;                
            x7_obs_new = x7_obs_new + Apow1 * delta_xN;
            x6_obs_new = x6_obs_new + Apow2 * delta_xN;
            x5_obs_new = x5_obs_new + Apow3 * delta_xN;
            x4_obs_new = x4_obs_new + Apow4 * delta_xN;
            x3_obs_new = x3_obs_new + Apow5 * delta_xN;
            x2_obs_new = x2_obs_new + Apow6 * delta_xN;
            x1_obs_new = x1_obs_new + Apow7 * delta_xN; 
            x0_obs_new = x0_obs_new + Apow8 * delta_xN; % -> result
        end
        
        xd_obs_old = xd_obs_new;
        x0_obs_old = x0_obs_new;
        x1_obs_old = x1_obs_new;
        x2_obs_old = x2_obs_new;
        x3_obs_old = x3_obs_new;
        x4_obs_old = x4_obs_new;
        x5_obs_old = x5_obs_new;
        x6_obs_old = x6_obs_new;
        x7_obs_old = x7_obs_new;
        if (n_delay==9); x8_obs_old = x8_obs_new; end
        
        obs_x0(:,k-n_delay) = x0_obs_new;
        obs_xd(:,k-n_delay) = xd_obs_new;

        % Compute q-vector
        fgm_x0(:,k-n_delay) = x0_obs_new;
        fgm_xd(:,k-n_delay) = xd_obs_new;
        fgm_u(:,k-n_delay) = single(u_sim(id_to_cm,k-1));
        q = q_mat*[x0_obs_new; xd_obs_new];

        % Compute lower and upper limit
        if use_single == true
            lower_u = max(-u_max-single(SOFB_setp),-u_rate+single(y_awr));
            upper_u = min(u_max-single(SOFB_setp),u_rate+single(y_awr));
        else
            lower_u = max(-u_max-double(SOFB_setp),-u_rate+double(y_awr));
            upper_u = min(u_max-double(SOFB_setp),u_rate+double(y_awr));
        end
        assert(sum(lower_u > upper_u)==0)

        % Fast gradient method
        if use_single == true
            out_global = single(u_sim(id_to_cm,k-1));
        else            
            out_global = double(u_sim(id_to_cm,k-1));
        end
        for i_iter = 1 : 1 : MAX_ITER
            z_old = z_new;
            t = J*out_global - q;
            z_new = max(lower_u, min(upper_u, t));
            out_global = (1+beta_fgm) * z_new - beta_fgm*z_old;
        end
        fgm_result = z_new; % new: use z_new instead of out_global
        fgm_out(:,k-n_delay) = fgm_result;
        if hil_mode == true
            fgm_result = round(fgm_result*1e3,0)/1e3; % note that MPC setpoint is in mA and we are streaming microA
        end
        u_sim(id_to_cm,k) = double(fgm_result);
        
        % AWR
        x_awr_old = x_awr_new;
        x_awr_new = A_awr*x_awr_old + B_awr*(double(fgm_result));
        y_awr = C_awr*x_awr_new + D_awr*(double(fgm_result));
    end
    
    % Plant
    x_sim_old = x_sim_new;
    x_sim_new = Ap*x_sim_old + Bp*(u_sim(:,k)); % note that the plant model accepts mA, but in D-I we need ending up with A
    x_sim(:, k) = x_sim_old;
end
u_sim = u_sim';
y_sim = y_sim';
x_sim = x_sim';

end