%sim_mpc_OSQP.m Run MPC simulation using OSQP
%
% Inputs:
%   n_samples   : Number of samples. Must match size(dist,2)
%   n_delay     : Number of delay steps. Only 8 or 9 supported.
%   dist        : Disturbance
%   Ap, Bp, Cp  : Plant
%   Ao, Bo, Co  : State observer
%   Ad, Cd      : Disturbance observer
%   LxN_obs     : State observer gain
%   Lxd_obs     : Disturbance observer gain
%   J_MPC       : Hessian of QP 
%   q_mat       : Matrix to compute q = q_mat*[x0_obs_new; xd_obs_new]
%   y_max       : Output constraint amplitude limit
%   u_max       : Input constraint amplitude limit
%   u_rate      : Input constraint rate limit
%   id_to_bpm   : Controlled outputs selected as y(id_to_bpm)
%   id_to_cm    : Controlled inputs selected as u(id_to_cm)
%   A_awr, B_awr, C_awr, D_awr : State-space system for rate computation
%   SOFB_setp   : Existing setpoint for amplitude computation
%   ol_mode     : Run in open-loop mode
%
% Outputs:
%   y_sim       : Simulated output
%   u_sim       : Simulated inputs
%
function [y_sim,u_sim] = sim_mpc_OSQP(...
            n_samples, n_delay, dist,...
            Ap, Bp, Cp,... % Plant
            Ao, Bo, Co, Ad, Cd, LxN_obs, Lxd_obs,... % Observer
            J_MPC, q_mat, y_max,...
            u_max, u_rate,... % Input constraints
            id_to_bpm, id_to_cm,... % Mapping for used BPMs & CMs
            A_awr, B_awr, C_awr, D_awr,... % Used for rate constraints
            SOFB_setp,... % Additional magnet setpoints
            ol_mode)  % Run in open-loop mode
        
if ~exist('ol_mode','var')
    ol_mode = false;
end
        
%%
assert((n_delay == 8)||n_delay==9);
use_single = false;
hil_mode = true;


[nx_plant, nu_plant] = size(Bp);
[ny_plant, ~] = size(Cp);
[nx_obs, nu_obs] = size(Bo);
[ny_obs, ~] = size(Co);

% Variables Plant
x_sim_new=zeros(nx_plant,1); x_sim_old=zeros(nx_plant,1);
y_sim = zeros(ny_plant, n_samples);
u_sim = zeros(nu_plant, n_samples);

% Variables AWR
[ny_awr,nx_awr]=size(C_awr);
x_awr_new=zeros(nx_awr,1);
y_awr=zeros(ny_awr,1);

% Variables Observer
if use_single == true
    x0_obs_new=single(zeros(nx_obs,1));
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
else
    x0_obs_new=double(zeros(nx_obs,1));
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

% Variables Solver
if use_single == true
    J     = single(J_MPC);
    q_mat = single(q_mat);
    y_max = single(y_max);
    u_max = single(u_max);
    y_mat = single(Co*Ao);
else
    J     = double(J_MPC);
    q_mat = double(q_mat);
    y_max = double(y_max);
    u_max = double(u_max);
    y_mat = double(Co*Ao);
end


MAX_ITER = 20;
osqp_solver = osqp;
settings = osqp_solver.default_settings();
settings.verbose = 0;
settings.polish = 0;
settings.adaptive_rho = 0;
settings.max_iter = MAX_ITER;
settings.check_termination = MAX_ITER;

A_constr = [eye(ny_obs); Co*Bo];
l_constr = [max(-u_max-SOFB_setp,-u_rate+y_awr);...
            -y_max-y_mat*x0_obs_new];
u_constr = [min(u_max-SOFB_setp,u_rate+y_awr);...
            y_max-y_mat*x0_obs_new];
osqp_solver.setup(J, zeros(1,size(q_mat,1)), A_constr, l_constr, u_constr, settings);

x0_obs = zeros(nx_obs,n_samples);
xd_obs = zeros(ny_obs,n_samples);
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
        % Observer - state update
        if use_single == true
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
        
        x0_obs(:,k) = x0_obs_new;
        xd_obs(:,k) = xd_obs_new;
        
        % Compute q-vector
        q = q_mat*[x0_obs_new; xd_obs_new];

        % Compute lower and upper limit
        if use_single == true
            l_constr = [max(-u_max-single(SOFB_setp),-u_rate+single(y_awr));...
                        -y_max-y_mat*single(x0_obs_new)];
            u_constr = [min(u_max-SOFB_setp,u_rate+single(y_awr));...
                        y_max-y_mat*single(x0_obs_new)];
        else
            l_constr = [max(-u_max-double(SOFB_setp),-u_rate+double(y_awr));...
                        -y_max-y_mat*double(x0_obs_new)-double(xd_obs_new)];
            u_constr = [min(u_max-SOFB_setp,u_rate+double(y_awr));...
                        y_max-y_mat*double(x0_obs_new)-double(xd_obs_new)];
        end        
        osqp_solver.update('q', q, 'l', l_constr, 'u', u_constr);
        result = osqp_solver.solve();
        osqp_result = result.x(1:nu_obs);
        assert((result.info.status_val == 1) || (result.info.status_val == 2));
        if hil_mode == true
            osqp_result = round(osqp_result*1e6,0)/1e6;
        end
        u_sim(id_to_cm,k) = double(osqp_result);
        
        % AWR
        x_awr_old = x_awr_new;
        x_awr_new = A_awr*x_awr_old + B_awr*(double(osqp_result));
        y_awr = C_awr*x_awr_new + D_awr*(double(osqp_result));
    end
    
    % Plant
    x_sim_old = x_sim_new;
    x_sim_new = Ap*x_sim_old + Bp*(u_sim(:,k));
end
u_sim = u_sim';
y_sim = y_sim';

end