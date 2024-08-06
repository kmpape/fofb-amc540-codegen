addpath('/home/idris/Documents/EngSci/Matlab/my_functions');
addpath('/home/idris/Documents/EngSci/Matlab/models');
addpath('/home/idris/Documents/EngSci/Matlab/simulink_models');
addpath('/home/idris/Documents/EngSci/Matlab/osqp/osqp-0.4.1-matlab-linux64');
addpath('/home/idris/Documents/EngSci/Matlab');
addpath('..')
clear all
close all
clc

%% Options
print_csv = false;
fname_RM = '../ORMS/GoldenBPMResp_DIAD.mat';
fname_X = '../../DATA/04092022_135252_data_X.mat';
fname_Y = '../../DATA/04092022_135252_data_Y.mat';
pick_dir = 2;
dirs = {'horizontal','vertical'};
pick_direction = dirs{pick_dir};
do_step = false;
sim_IMC = false;
use_FGM = true;

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
fname = sprintf('mpc_data_13092022_nd%d.mat',n_delay);
if true % ~exist(fname,'file')
    print_msg = false;
    [Ao_x, Bo_x, Co_x, Ap_x, Bp_x, Cp_x, Ad_x, Cd_x,...
          Kfd_x, Kfx_x, Kcx_x, Kcd_x, P_x, Rlqr_x, Qlqr_x,...
          Ao_y, Bo_y, Co_y, Ap_y, Bp_y, Cp_y, Ad_y, Cd_y,...
          Kfd_y, Kfx_y, Kcx_y, Kcd_y, P_y, Rlqr_y, Qlqr_y] =...
          observer_regulator(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,n_delay,fname,print_msg);
else
    load(fname);
end

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
Lx8_obs = Kfx;
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
q_mat = [q_mat_x0, q_mat_xd];

if use_FGM == true
    eigmax = max(eig(J_mpc)); 
    eigmin = min(eig(J_mpc));
    J_mpc = eye(size(J_mpc))-J_mpc/eigmax;
    beta_fgm = (sqrt(eigmax) - sqrt(eigmin)) / (sqrt(eigmax) + sqrt(eigmin));
    q_mat = [q_mat_x0, q_mat_xd]/eigmax;
    fprintf("GOOD: cond(J_mpc)=%e\n",cond(J_mpc));
end

%% Rate limiter on VME processors
a_awr = 2*2*pi;
g_awr_z = tf([a_awr/(a_awr+2/Ts),a_awr/(a_awr+2/Ts)],...
    [1,(a_awr-2/Ts)/(a_awr+2/Ts)],Ts,'Variable','z^-1');
sys_awr = eye(nu) * g_awr_z;

%% Measurement Data
imode = 1;
n_samples = 300+n_delay+1;
if strcmp(pick_direction, 'vertical')
    [UR,SR,VR] = svd(RMy);
else
    [UR,SR,VR] = svd(RMx);
end
mag_u = 10;

% Simulation
endt = (n_samples*Ts)-Ts;
Lsim = n_samples*Ts;
t= 0:Ts:endt;
SOFB_setp = zeros(nu,1);
ss_awr = ss(sys_awr);
open_loop = true;
idist_vec = [1,50,100];

inds = 1:1:length(t)-n_delay;
dist_bpm_1 = zeros(length(inds),length(idist_vec)+1);
estim_bpm_1 = zeros(length(inds),length(idist_vec)+1);
fig = figure; hold on;
for idist = 1:length(idist_vec)+1
    if idist > length(idist_vec)
        if strcmp(pick_direction, 'vertical')
            fname = fname_Y;
        else
            fname = fname_X;
        end
        distdata = load(fname);
        doff = -distdata.y(:,1:n_samples);
    else
        imode = idist_vec(idist);
        tmp = UR(:,imode)*mag_u/UR(1,imode); % *SR(imode,imode)
        doff_tmp = zeros(TOT_BPM,1);
        doff_tmp(id_to_bpm) = tmp;
        doff = doff_tmp .* ones(1,n_samples);
        doff = doff+randn(size(doff)); 
    end
    [y_sim,u_sim,~,~,~,obs_x0,obs_xd,fgm_x0,fgm_xd,fgm_u,fgm_out,~,~] = sim_mpc(...
                    n_samples, n_delay, doff,...
                    Ap, Bp, Cp, ... % Plant
                    Ao, Bo, Co, Ad, Cd, Lx8_obs, Lxd_obs,... % Observer
                    J_mpc , q_mat, beta_fgm,... % FGM
                    u_max , u_rate,... % FGM
                    id_to_bpm, id_to_cm,...
                    ss_awr.A,ss_awr.B,ss_awr.C,ss_awr.D,...
                    SOFB_setp,open_loop);
    scale_u = 1e-3;
    if false
        subplot(length(idist_vec)+1,2,(idist-1)*2+1); plot(t(inds),doff(id_to_bpm,inds)'); title('Disturbance MPC');
        subplot(length(idist_vec)+1,2,(idist-1)*2+2); plot(t(inds),obs_xd(:,inds)'); title('Disturbance Estimate MPC');
    else
        subplot(length(idist_vec)+1,1,idist); 
        plot(t(inds),doff(1,inds)'); title('Disturbance MPC'); hold on;
        plot(t(inds),obs_xd(1,inds)'); title('Disturbance Estimate MPC');
    end
    dist_bpm_1(:,idist) = doff(1,inds);
    estim_bpm_1(:,idist) = obs_xd(1,inds);
end
set(fig, 'position',[100 100 1900 1600]);
if print_csv
    tt = t(inds);
    cols = {'x','distmode1','distmode50','distmode100','dist','distmode1hat','distmode50hat','distmode100hat','disthat'};
    data = [tt'*1000,dist_bpm_1,estim_bpm_1];
    csvwrite_with_headers('thesis_csv/obsstep.csv',data,cols);
end
