% Change path to OSQP here
addpath('/home/idris/Documents/EngSci/Matlab/osqp/osqp-0.4.1-matlab-linux64');
addpath('..')

%% Options
fname_RM = '../ORMS/GoldenBPMResp_DIAD.mat';
fname_X = '../DATA/04092022_135252_data_X.mat';
fname_Y = '../DATA/04092022_135252_data_Y.mat';
pick_dir = 2;
dirs = {'horizontal','vertical'};
pick_direction = dirs{pick_dir};
sim_IMC = false;
use_FGM = true;  % If false, use OSQP

%% Hardlimits
load('../ORMS/correctors.mat');
hardlimits = corrector_data.MaxAmps(1:172); % in Amperes

%% Configure Diamond-I Storage Ring
load(fname_RM);
RMorigx = Rmat(1).Data(:,:);%  * 1e6; % RM(1) eq to RM(1,1)
[ny_x, nu_x] = size(RMorigx);
RMorigy = Rmat(4).Data(:,:);%  * 1e6; % RM(4) eq to RM(2,2)
[ny_y, nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
[TOT_BPM, TOT_CM] = size(RMorigx);
square_config = true;
[id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v5(RMorigx,RMorigy,square_config);
RMx = RMorigx(id_to_bpm_x,id_to_cm_x);
RMy = RMorigy(id_to_bpm_y,id_to_cm_y);

%% Observer and Regulator
n_delay = 9;
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time [s]
% Obtain LQR and Observer matrices
[Ao_x, Bo_x, Co_x, Ap_x, Bp_x, Cp_x, Ad_x, Cd_x,...
      Kfd_x, Kfx_x, Kcx_x, Kcd_x, P_x, Rlqr_x, Qlqr_x,...
      Ao_y, Bo_y, Co_y, Ap_y, Bp_y, Cp_y, Ad_y, Cd_y,...
      Kfd_y, Kfx_y, Kcx_y, Kcd_y, P_y, Rlqr_y, Qlqr_y] =...
      observer_regulator(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,n_delay);

if strcmp(pick_direction, 'vertical')
    id_to_bpm = id_to_bpm_y;
    id_to_cm = id_to_cm_y;
    RM = RMy;
    aI_Hz = 700; % Corrector bandwidth [Hz]
    Ao = Ao_y; Bo = Bo_y; Co = Co_y; Ad = Ad_y; Cd = Cd_y; % plant for observer
    Ap = Ap_y; Bp = Bp_y; Cp = Cp_y; % plant with all BPMs and CMs
    Kfd = Kfd_y; % Observer gain for disturbance
    Kfx = Kfx_y; % Observer gain for state
    P_mpc = P_y; % terminal cost matrix
    Q_mpc = Qlqr_y; % state weighting matrix
    R_mpc = Rlqr_y; % input weighting matrix
else
    id_to_bpm = id_to_bpm_x;
    id_to_cm = id_to_cm_x;
    RM = RMx;
    aI_Hz = 500; % Corrector bandwidth [Hz]
    Ao = Ao_x; Bo = Bo_x; Co = Co_x; Ad = Ad_x; Cd = Cd_x; % plant for observer
    Ap = Ap_x; Bp = Bp_x; Cp = Cp_x; % plant with all BPMs and CMs
    Kfd = Kfd_x; % Observer gain for disturbance
    Kfx = Kfx_x; % Observer gain for state
    P_mpc = P_x; % terminal cost matrix
    Q_mpc = Qlqr_x; % state weighting matrix
    R_mpc = Rlqr_x; % input weighting matrix
end
[ny, nu] = size(RM);
nx = nu;


%% Observer
Lxd_obs = Kfd;
Lx8_obs = Kfx;  % NOTE: this should be zero
S_sp_pinv = pinv([eye(nx)-Ao, -Bo; Co, zeros(ny, nu)]);
S_sp_pinv = S_sp_pinv(:, nx+1:end);

%% MPC
horizon = 1;
u_rate_scalar = 1*1000;
u_rate = u_rate_scalar*ones(nu,1);
u_max = hardlimits(id_to_cm)*1000;
J_mpc = Bo'*P_mpc*Bo+R_mpc;
S_sp_pinv_x = S_sp_pinv(1:nx,:);
S_sp_pinv_u = S_sp_pinv(nx+1:end,:);
q_mat_x0 = Bo'*P_mpc*Ao;
q_mat_xd = -[Bo'*P_mpc, R_mpc]*[S_sp_pinv_x; S_sp_pinv_u]*(-Cd);
q_mat = [q_mat_x0, q_mat_xd];  % Matrix to obain q(x_obs)

if use_FGM == true
    eigmax = max(eig(J_mpc)); 
    eigmin = min(eig(J_mpc));
    J_mpc = eye(size(J_mpc))-J_mpc/eigmax;
    beta_fgm = (sqrt(eigmax) - sqrt(eigmin)) / (sqrt(eigmax) + sqrt(eigmin));
    q_mat = [q_mat_x0, q_mat_xd]/eigmax;
end

%% Rate limiter on VME processors
a_awr = 2*2*pi;
g_awr_z = tf([a_awr/(a_awr+2/Ts),a_awr/(a_awr+2/Ts)],...
    [1,(a_awr-2/Ts)/(a_awr+2/Ts)],Ts,'Variable','z^-1');
sys_awr = eye(nu) * g_awr_z;
ss_awr = ss(sys_awr);

%% Measurement Data
imode = 1; % select mode for disturbance
n_samples = 100;
if strcmp(pick_direction, 'vertical')
    [UR,SR,VR] = svd(RMy);
else
    [UR,SR,VR] = svd(RMx);
end
mag_u = 10;
tmp = UR(:,imode)*SR(imode,imode)*mag_u;
doff_tmp = zeros(TOT_BPM,1);
doff_tmp(id_to_bpm) = tmp;
doff = doff_tmp .* ones(1,n_samples);

% Simulation
endt = (n_samples*Ts)-Ts;
Lsim = n_samples*Ts;
t= 0:Ts:endt;

SOFB_setp = zeros(nu,1);

%% Run MPC
if use_FGM == true
    ind = find(id_to_cm==38);
    [y_sim,u_sim,~,obs_y,obs_u,obs_x0,obs_xd,fgm_x0,fgm_xd,fgm_u,fgm_out,lower_u,upper_u] = sim_mpc(...
                    n_samples, n_delay, doff,...
                    Ap, Bp, Cp, ... % Plant
                    Ao, Bo, Co, Ad, Cd, Lx8_obs, Lxd_obs,... % Observer
                    J_mpc , q_mat, beta_fgm,... % FGM
                    u_max , u_rate,... % FGM
                    id_to_bpm, id_to_cm,...
                    ss_awr.A,ss_awr.B,ss_awr.C,ss_awr.D,...
                    SOFB_setp,false);
else
    y_max_scalar = Inf;
    y_max = ones(length(id_to_bpm),1)*y_max_scalar;
    [y_sim,u_sim] = sim_mpc_OSQP(...
                n_samples, n_delay, doff,...
                Ap, Bp, Cp, ... % Plant
                Ao, Bo, Co, Ad, Cd, Lx8_obs, Lxd_obs,... % Observer
                J_mpc , q_mat, y_max,... % FGM
                u_max , u_rate,... % FGM
                id_to_bpm, id_to_cm,...
                ss_awr.A,ss_awr.B,ss_awr.C,ss_awr.D,...
                SOFB_setp,false);
end

%% Run IMC
if sim_IMC
    addpath('../IMC')
    
    bw = 1/(n_delay*Ts);
    [network_scaling,...
      Acx, Bcx, Ccx, Dcx, Ax, Bx, Cx, Dx,...
      Acy, Bcy, Ccy, Dcy, Ay, By, Cy, Dy] = get_IMC_controller_cfgv4_square(RMorigx, RMorigy,...
                                                                                bw, n_delay);
    if strcmp(pick_direction, 'vertical')
        Ac_imc = Acy; Bc_imc = Bcy; Cc_imc = Ccy; Dc_imc = Dcy;
        Ap_imc = Ay; Bp_imc = By; Cp_imc = Cy; Dp_imc = Dy;
    else        
        Ac_imc = Acx; Bc_imc = Bcx; Cc_imc = Ccx; Dc_imc = Dcx;
        Ap_imc = Ax; Bp_imc = Bx; Cp_imc = Cx; Dp_imc = Dx;
    end
    [y_sim_imc, u_sim_imc] = sim_standard_imc(...
            n_samples, n_delay, id_to_bpm, id_to_cm, doff,...
            Ac_imc, Bc_imc, Cc_imc, Dc_imc,... % CONTROLLER STATE-SPACE
            Ap_imc, Bp_imc, Cp_imc, Dp_imc, network_scaling); % PLANT STATE-SPACE                                                                            
end

scale_u = 1e-3;

y_awr = lsim(sys_awr,u_sim(1:length(t),id_to_cm)*scale_u,t);
if sim_IMC
    figure;
    subplot(2,2,1); plot(doff(id_to_bpm,:)'); title('Disturbance IMC');
    subplot(2,2,2); plot(y_sim_imc(:,id_to_bpm)); title('Output IMC');
    subplot(2,2,3); plot(u_sim_imc*scale_u); title('Input IMC');
end

%% Run LQR
if strcmp(pick_direction, 'vertical')
    id_to_bpm = id_to_bpm_y;
    id_to_cm = id_to_cm_y;
    RM = RMy;
    RMorig = RMorigy;
    aI_Hz = 700; % Corrector bandwidth [Hz]
    Ao = Ao_y; Bo = Bo_y; Co = Co_y; Ad = Ad_y; Cd = Cd_y; % plant for observer
    Ap = Ap_y; Bp = Bp_y; Cp = Cp_y; % plant with all BPMs and CMs
    Kfd = Kfd_y; % Observer gain for disturbance
    Kfx = Kfx_y; % Observer gain for state
    Kcd = Kcd_y;
    Kcx = Kcx_y;
else
    id_to_bpm = id_to_bpm_x;
    id_to_cm = id_to_cm_x;
    RM = RMx;
    RMorig = RMorigx;
    aI_Hz = 500; % Corrector bandwidth [Hz]
    Ao = Ao_x; Bo = Bo_x; Co = Co_x; Ad = Ad_x; Cd = Cd_x; % plant for observer
    Ap = Ap_x; Bp = Bp_x; Cp = Cp_x; % plant with all BPMs and CMs
    Kfd = Kfd_x; % Observer gain for disturbance
    Kfx = Kfx_x; % Observer gain for state
    Kcd = Kcd_x;
    Kcx = Kcx_x;
end

Astate = [Ao, zeros(nu,nu*n_delay); eye(nu*n_delay), zeros(nu*n_delay,nu)];
Bstate = [Bo; zeros(nu*n_delay,nu)];
Cstate = [zeros(ny,nu*n_delay), Co];
Kfstate = [zeros(nu*n_delay, ny); Kfx];
for i=1:n_delay; Kfstate(1+(i-1)*nu:i*nu,:) = Ao^(n_delay+1-i)*Kfx; end
Kcstate = [Kcx, zeros(nu,nu*n_delay)];

Aobs = blkdiag(Astate, Ad);
Bobs = [Bstate; zeros(ny,nu)];
Cobs = [Cstate, Cd];
Kc = [Kcstate, Kcd];
Kf = [Kfstate; Kfd];
[y_sim_lqr, u_sim_lqr, ~] = sim_lqr_w_constraints(...
            n_samples, n_delay, doff,...
            Ap, Bp, Cp,... % plant
            Aobs, Bobs, Cobs, Kf, Kc,... % observer and regulator
            id_to_bpm, id_to_cm,...
            u_max , u_rate,... % FGM
            ss_awr.A,ss_awr.B,ss_awr.C,ss_awr.D,...
            false);
scale_u = 1e-3;
y_awr_lqr = lsim(sys_awr,u_sim_lqr(1:length(t),id_to_cm)*scale_u,t);

fig = figure;
subplot(2,4,1); plot(doff(id_to_bpm,:)'); title('Disturbance MPC');
subplot(2,4,2); plot(y_sim(:,id_to_bpm)); title('Output MPC');
subplot(2,4,3); plot(u_sim*scale_u); title('Input MPC');
subplot(2,4,4); plot(u_sim(1:length(t),id_to_cm)*scale_u-y_awr); title('AWR MPC');
subplot(2,4,5); plot(doff(id_to_bpm,:)'); title('Disturbance LQR');
subplot(2,4,6); plot(y_sim_lqr(:,id_to_bpm)); title('Output LQR');
subplot(2,4,7); plot(u_sim_lqr*scale_u); title('Input LQR');
subplot(2,4,8); plot(u_sim_lqr(1:length(t),id_to_cm)*scale_u-y_awr_lqr); title('AWR LQR');
set(fig, 'position',[100 100 1900 1600]);

doff_mode = UR'*doff(id_to_bpm,:);
ysim_mode = y_sim(:,id_to_bpm)*UR;
u_sim_mode = u_sim(:,id_to_cm)*VR*scale_u;
ysim_lqr_mode = y_sim_lqr(:,id_to_bpm)*UR;
u_sim_lqr_mode = u_sim_lqr(:,id_to_cm)*VR*scale_u;
fig = figure;
subplot(2,4,1); plot(doff_mode'); title('MODE Disturbance MPC');
subplot(2,4,2); plot(ysim_mode); title('MODE Output MPC');
subplot(2,4,3); plot(u_sim_mode); title('MODE Input MPC');
subplot(2,4,4); plot(u_sim(1:length(t),id_to_cm)*scale_u-y_awr); title('AWR MPC');
subplot(2,4,5); plot(doff_mode'); title('MODE Disturbance LQR');
subplot(2,4,6); plot(ysim_lqr_mode); title('MODE Output LQR');
subplot(2,4,7); plot(u_sim_lqr_mode); title('MODE Input LQR');
set(fig, 'position',[100 100 1900 1600]);

scale_mode = 1e-3;
