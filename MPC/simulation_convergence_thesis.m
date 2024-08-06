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
if strcmp(pick_direction, 'vertical')
    fname_dist = '../../DATA/04092022_135000_data_Y.mat';
else
    fname_dist = '../../DATA/04092022_135000_data_X.mat';
end
distdata = load(fname_dist);
n_samples = 10000+n_delay+1;
%
Rsettings = {'a','b','c','d'};
n_iters_all = zeros(n_samples,length(Rsettings));
abs_err_qp_all = zeros(n_samples,length(Rsettings));
rel_err_qp_all = zeros(n_samples,length(Rsettings));
for isetting = 1:length(Rsettings)
    Rsetting = Rsettings{isetting};
    print_msg = false;
    [Ao_x, Bo_x, Co_x, Ap_x, Bp_x, Cp_x, Ad_x, Cd_x,...
          Kfd_x, Kfx_x, Kcx_x, Kcd_x, P_x, Rlqr_x, Qlqr_x,...
          Ao_y, Bo_y, Co_y, Ap_y, Bp_y, Cp_y, Ad_y, Cd_y,...
          Kfd_y, Kfx_y, Kcx_y, Kcd_y, P_y, Rlqr_y, Qlqr_y] =...
          observer_regulator_bad(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,n_delay,'',print_msg,Rsetting);
end
asdf
for i=1:1
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

    % Observer
    Lxd_obs = Kfd;
    Lx8_obs = Kfx;
    S_sp_pinv = pinv([eye(nx)-Ao, -Bo; Co, zeros(ny, nu)]);
    S_sp_pinv = S_sp_pinv(:, nx+1:end);

    % MPC
    u_rate_scalar = 1*1000;
    u_rate = u_rate_scalar*ones(nu,1);
    u_max = hardlimits(id_to_cm)*1000;
    J_mpc = Bo'*P_mpc*Bo+R_mpc; J_mpc = 0.5*(J_mpc+J_mpc');
    fprintf("SETTING=%s: cond(R_mpc)=%e, cond(P_mpc)=%e, cond(J_mpc)=%e\n",Rsetting,cond(R_mpc),cond(P_mpc),cond(J_mpc));
    JJ = J_mpc;
    S_sp_pinv_x = S_sp_pinv(1:nx,:);
    S_sp_pinv_u = S_sp_pinv(nx+1:end,:);
    q_mat_x0 = Bo'*P_mpc*Ao;
    q_mat_xd = -[Bo'*P_mpc, R_mpc]*[S_sp_pinv_x; S_sp_pinv_u]*(-Cd);
    qq_mat = [q_mat_x0, q_mat_xd];
    eigmax = max(eig(J_mpc)); 
    eigmin = min(eig(J_mpc));
    J_mpc = eye(size(J_mpc))-J_mpc/eigmax;
    beta_fgm = (sqrt(eigmax) - sqrt(eigmin)) / (sqrt(eigmax) + sqrt(eigmin));
    q_mat = [q_mat_x0, q_mat_xd]/eigmax;

    % Rate limiter on VME processors
    a_awr = 2*2*pi;
    g_awr_z = tf([a_awr/(a_awr+2/Ts),a_awr/(a_awr+2/Ts)],...
        [1,(a_awr-2/Ts)/(a_awr+2/Ts)],Ts,'Variable','z^-1');
    sys_awr = eye(nu) * g_awr_z;

    % Measurement Data
    doff = distdata.y(:,1:n_samples);

    % Simulation
    endt = (n_samples*Ts)-Ts;
    Lsim = n_samples*Ts;
    t= 0:Ts:endt;
    SOFB_setp = zeros(nu,1);
    ss_awr = ss(sys_awr);
    open_loop = false;
    [y_sim,u_sim,x_sim,n_iters,abs_err_qp,rel_err_qp] = sim_mpc_convergence(...
                        n_samples, n_delay, doff,...
                        Ap, Bp, Cp, ... % Plant
                        Ao, Bo, Co, Ad, Cd, Lx8_obs, Lxd_obs,... % Observer
                        J_mpc , q_mat, beta_fgm,... % FGM
                        u_max , u_rate,... % FGM
                        id_to_bpm, id_to_cm,...
                        ss_awr.A,ss_awr.B,ss_awr.C,ss_awr.D,...
                        SOFB_setp,(JJ+JJ')/2,qq_mat);
    n_iters_all(:,isetting) = n_iters';
    abs_err_qp_all(:,isetting) = abs_err_qp';
    rel_err_qp_all(:,isetting) = rel_err_qp';
    
    
end
for isetting = 1:length(Rsettings)
    Rsetting = Rsettings{isetting};
	figure;
    subplot(3,1,1)
    plot(n_iters_all(:,isetting))
    title(Rsetting)
    subplot(3,1,2)
    plot(abs_err_qp_all(:,isetting))
    subplot(3,1,3)
    plot(rel_err_qp_all(:,isetting))
end
if true
    indsp = 100:1:size(n_iters_all,1)-10;
    it_mean = mean(n_iters_all(indsp,:),1);
    it_std = std(n_iters_all(indsp,:));
%     it_data = [it_mean',it_mean'-it_std',it_mean'+it_std',...
%         min(n_iters_all(indsp,:),[],1)',max(n_iters_all(indsp,:),[],1)',...
%         median(n_iters_all(indsp,:),1)'];
    it_data = [it_mean',it_std',it_std',...
        min(n_iters_all(indsp,:),[],1)',max(n_iters_all(indsp,:),[],1)',...
        median(n_iters_all(indsp,:),1)'];
    abs_mean = mean(abs_err_qp_all(indsp,:),1);
    abs_std = std(abs_err_qp_all(indsp,:));
%     abs_data = [abs_mean',abs_mean'-abs_std',abs_mean'+abs_std',...
%         min(abs_err_qp_all(indsp,:),[],1)',max(abs_err_qp_all(indsp,:),[],1)',...
%         median(abs_err_qp_all(indsp,:),1)'];
    abs_data = [abs_mean',abs_std',abs_std',...
        min(abs_err_qp_all(indsp,:),[],1)',max(abs_err_qp_all(indsp,:),[],1)',...
        median(abs_err_qp_all(indsp,:),1)'];
    rel_mean = mean(rel_err_qp_all(indsp,:),1);
    rel_std = std(rel_err_qp_all(indsp,:));
    rel_data = [rel_mean',rel_std',rel_std',...
        min(rel_err_qp_all(indsp,:),[],1)',max(rel_err_qp_all(indsp,:),[],1)',...
        median(rel_err_qp_all(indsp,:),1)'];
    
     
    x = [1,2,3]';
    cols = {'x',...
                'it_mean','it_mean_m','it_mean_p','it_min','it_max','it_median',...
                'abs_mean','abs_mean_m','abs_mean_p','abs_min','abs_max','abs_median',...
                'rel_mean','rel_mean_m','rel_mean_p','rel_min','rel_max','rel_median'};
    data = [x,it_data,abs_data,rel_data*100];
    csvwrite_with_headers('thesis_csv/convergence_20_eps1e-3orInf.csv',data,cols);
end


%% Performance
n_delay = 9;
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time [s]
if strcmp(pick_direction, 'vertical')
    fname_dist = '../../DATA/04092022_135000_data_Y.mat';
else
    fname_dist = '../../DATA/04092022_135000_data_X.mat';
end
distdata = load(fname_dist);
n_samples = 100700;
%
Rsettings = {'b','c','d'};
for isetting = 1:length(Rsettings)
    Rsetting = Rsettings{isetting};
    print_msg = false;
    [Ao_x, Bo_x, Co_x, Ap_x, Bp_x, Cp_x, Ad_x, Cd_x,...
          Kfd_x, Kfx_x, Kcx_x, Kcd_x, P_x, Rlqr_x, Qlqr_x,...
          Ao_y, Bo_y, Co_y, Ap_y, Bp_y, Cp_y, Ad_y, Cd_y,...
          Kfd_y, Kfx_y, Kcx_y, Kcd_y, P_y, Rlqr_y, Qlqr_y] =...
          observer_regulator_bad(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,n_delay,'',print_msg,Rsetting);
    
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

    % Observer
    Lxd_obs = Kfd;
    Lx8_obs = Kfx;
    S_sp_pinv = pinv([eye(nx)-Ao, -Bo; Co, zeros(ny, nu)]);
    S_sp_pinv = S_sp_pinv(:, nx+1:end);

    % MPC
    u_rate_scalar = 1*1000;
    u_rate = u_rate_scalar*ones(nu,1);
    u_max = hardlimits(id_to_cm)*1000;
    J_mpc = Bo'*P_mpc*Bo+R_mpc; J_mpc = 0.5*(J_mpc+J_mpc');
    fprintf("SETTING=%s: cond(R_mpc)=%e, cond(P_mpc)=%e, cond(J_mpc)=%e\n",Rsetting,cond(R_mpc),cond(P_mpc),cond(J_mpc));
    JJ = J_mpc;
    S_sp_pinv_x = S_sp_pinv(1:nx,:);
    S_sp_pinv_u = S_sp_pinv(nx+1:end,:);
    q_mat_x0 = Bo'*P_mpc*Ao;
    q_mat_xd = -[Bo'*P_mpc, R_mpc]*[S_sp_pinv_x; S_sp_pinv_u]*(-Cd);
    qq_mat = [q_mat_x0, q_mat_xd];
    eigmax = max(eig(J_mpc)); 
    eigmin = min(eig(J_mpc));
    J_mpc = eye(size(J_mpc))-J_mpc/eigmax;
    beta_fgm = (sqrt(eigmax) - sqrt(eigmin)) / (sqrt(eigmax) + sqrt(eigmin));
    q_mat = [q_mat_x0, q_mat_xd]/eigmax;

    % Rate limiter on VME processors
    a_awr = 2*2*pi;
    g_awr_z = tf([a_awr/(a_awr+2/Ts),a_awr/(a_awr+2/Ts)],...
        [1,(a_awr-2/Ts)/(a_awr+2/Ts)],Ts,'Variable','z^-1');
    sys_awr = eye(nu) * g_awr_z;

    % Measurement Data
    doff = distdata.y(:,1:n_samples);

    % Simulation
    endt = (n_samples*Ts)-Ts;
    Lsim = n_samples*Ts;
    t= 0:Ts:endt;
    SOFB_setp = zeros(nu,1);
    ss_awr = ss(sys_awr);
    open_loop = false;
    [y_sim,u_sim,~] = sim_mpc_convergence_fast(...
                        n_samples, n_delay, doff,...
                        Ap, Bp, Cp, ... % Plant
                        Ao, Bo, Co, Ad, Cd, Lx8_obs, Lxd_obs,... % Observer
                        J_mpc , q_mat, beta_fgm,... % FGM
                        u_max , u_rate,... % FGM
                        id_to_bpm, id_to_cm,...
                        ss_awr.A,ss_awr.B,ss_awr.C,ss_awr.D,...
                        SOFB_setp,(JJ+JJ')/2,qq_mat,false);
    y_sim = y_sim(500:end,:);
    nfft = 100000;
    Fs = 10072;
    [psdy,Fhz] = pwelch(y_sim,nfft/10,[],nfft,Fs);
    df = Fhz(2);
    w_psd = [(0.1:0.1:1), (1:1:10), (15:5:100), (110:10:1000), (1100:100:10000)]*df; w_psd(w_psd > 5000) = [];
    inds_psd = get_log_from_linear_inds(w_psd, Fhz);
    figure; loglog(Fhz(inds_psd),(sqrt(abs(psdy(inds_psd,1))))); hold on;
    loglog(Fhz(:),(sqrt(abs(psdy(:,1)))));
    
    cols = {'freq', 'yasd_1', 'yasd_2', 'yasd_3', 'yasd_4', 'yasd_5', 'yasd_6', 'yasd_7'};
    data = [Fhz(inds_psd),sqrt(abs(psdy(inds_psd,1:7)))];
    csvwrite_with_headers(sprintf('thesis_csv/psdy_20_eps1e-3orInf_%s.csv',Rsetting),data,cols);
end
