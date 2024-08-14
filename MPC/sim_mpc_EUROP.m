%sim_mpc.m Run MPC simulation using fast gradient method
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
%   J_MPC       : Modified Hessian of QP in FGM form
%   beta_fgm    : FGM step size
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
%   x_sim       : Simulated states
%   lower_u     : Lower bound for u
%   upper_u     : Lower bound for u
%
function [y_sim,u_sim,x_sim,lower_u,upper_u] = sim_mpc_EUROP(...
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
[nx_plant, nu_plant] = size(Bp);
[ny_plant, ~] = size(Cp);
[nx_obs, nu_obs] = size(Bo);
[ny_obs, ~] = size(Co);

% Variables Plant
x_sim_new = zeros(nx_plant,1); 
x_sim_old = zeros(nx_plant,1);
y_sim = zeros(ny_plant, n_samples);
u_sim = zeros(nu_plant, n_samples);
x_sim = zeros(nu_plant, n_samples);

% Variables AWR
[ny_awr,nx_awr] = size(C_awr);
x_awr_new = zeros(nx_awr,1);
y_awr = zeros(ny_awr,1);

% Variables Observer
x_obs_old = zeros(nx_obs, n_delay+1);
x_obs_new = zeros(nx_obs, n_delay+1);
xd_obs_old = zeros(ny_obs, 1);
ApowN = zeros(nx_obs*(n_delay+1), nx_obs);
for i = 1 : n_delay+1
    ApowN(1+(i-1)*nx_obs:i*nx_obs, :) = Ao.^(i-1);
end

% Variables Fast Gradient
J = J_MPC;
z_new=zeros(nu_obs,1);

MAX_ITER = 20;
for k = 1:1:n_samples

    % Measurement
    if ~ol_mode
        y_sim(:, k) = Cp*x_sim_new + dist(:, k);
    else
        y_sim(:, k) = dist(:, k);
    end

    if k > n_delay
        y_meas = y_sim(id_to_bpm, k-n_delay);
        x_obs_new(:,1) = Ao*x_obs_old(:,1) + Bo*u_sim(id_to_cm,k-1);
        x_obs_new(:,2:end) = x_obs_old(:,1:end-1);
        xd_obs_new = Ad*xd_obs_old;

        % Observer - measurement update
        delta_y = y_meas - Co*x_obs_new(:,end) - Cd*xd_obs_new;
        delta_xN = LxN_obs * delta_y;
        delta_xd = Lxd_obs * delta_y;
        xd_obs_new = xd_obs_new + delta_xd;
        x_obs_new = x_obs_new +fliplr(reshape(ApowN * delta_xN,...
                                              [nx_obs, n_delay+1]));
        
        %Copy new to old state
        xd_obs_old = xd_obs_new;
        x_obs_old(:,1:end-1) = x_obs_new(:,1:end-1);



        % Compute q vector
        q = q_mat*[x_obs_new(:, 1); xd_obs_new];

        % Compute lower and upper limit
        lower_u = max(-u_max-double(SOFB_setp),-u_rate+y_awr);
        upper_u = min(u_max-double(SOFB_setp),u_rate+y_awr);
        assert(sum(lower_u > upper_u)==0)

        % Fast gradient method
        out_global = u_sim(id_to_cm,k-1);
        for i_iter = 1 : 1 : MAX_ITER
            z_old = z_new;
            t = J*out_global - q;
            z_new = max(lower_u, min(upper_u, t));
            out_global = (1+beta_fgm) * z_new - beta_fgm*z_old;
        end
        fgm_result = z_new; % new: use z_new instead of out_global
        u_sim(id_to_cm,k) = fgm_result;
        
        % AWR
        x_awr_old = x_awr_new;
        x_awr_new = A_awr*x_awr_old + B_awr*fgm_result;
        y_awr = C_awr*x_awr_new + D_awr*fgm_result;
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