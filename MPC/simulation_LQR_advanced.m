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
fname_RM = '../ORMS/GoldenBPMResp_DIAD.mat';
fname_X = '../../DATA/04092022_135000_data_X.mat';
fname_Y = '../../DATA/04092022_135000_data_Y.mat';

pick_dir = 2;
dirs = {'horizontal','vertical'};
pick_direction = dirs{pick_dir};
do_step = false;
sim_IMC = false;

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
[id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v4(RMorigx,RMorigy,true);
RMx = RMorigx(id_to_bpm_x,id_to_cm_x);
RMy = RMorigy(id_to_bpm_y,id_to_cm_y);

%% Observer and Regulator
n_delay = 9;
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time [s]
fname = sprintf('mpc_data_02092022_nd%d.mat',n_delay);
if ~exist(fname,'file')
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
    Kcd = Kcd_y;
    Kcx = Kcx_y;
else
    id_to_bpm = id_to_bpm_x;
    id_to_cm = id_to_cm_x;
    RM = RMx;
    aI_Hz = 500; % Corrector bandwidth [Hz]
    Ao = Ao_x; Bo = Bo_x; Co = Co_x; Ad = Ad_x; Cd = Cd_x; % plant for observer
    Ap = Ap_x; Bp = Bp_x; Cp = Cp_x; % plant with all BPMs and CMs
    Kfd = Kfd_x; % Observer gain for disturbance
    Kfx = Kfx_x; % Observer gain for state
    Kcd = Kcd_x;
    Kcx = Kcx_x;
end
[ny, nu] = size(RM);
nx = nu;

%% Observer
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

%% Mode by mode observer
fname_mode = 'observer_advanced.mat';
if ~exist(fname_mode,'file')
    [Aoi_x, Boi_x, Coi_x, Kfi_x, Kci_x, Aoi_y, Boi_y, Coi_y, Kfi_y, Kci_y] =...
        observer_regulator_v2(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,n_delay,fname_mode,false);
else
    load(fname_mode);
end

if strcmp(pick_direction, 'vertical')
    Aoi = Aoi_y; Boi = Boi_y; Coi = Coi_y; Kfi = Kfi_y; Kci = Kci_y;
else
    Aoi = Aoi_x; Boi = Boi_x; Coi = Coi_x; Kfi = Kfi_x; Kci = Kci_x;
end
[UR,SR,VR] = svd(RM,'econ');

%% Measurement Data
if do_step
    n_samples = 2000;
    doff = ones(TOT_BPM,1) .* ones(1,n_samples)*10;
    don = doff;    
else
    if strcmp(pick_direction, 'vertical')
        fname = fname_Y;
    else
        fname = fname_X;
    end
    inds = 1:50000;
    load(fname);
    doff = y(:,inds);
    u_FOFB = u(:,inds);
    n_samples = length(inds);
    don = doff;
    all_bpms = 1:1:TOT_BPM;
    bad_bpm = setdiff(all_bpms,id_to_bpm);
    doff(bad_bpm,:) = 0;
    don(bad_bpm,:) = 0;
end

%% Simulation
endt = (n_samples*Ts)-Ts;
Lsim = n_samples*Ts;
t= 0:Ts:endt;

[y_sim, u_sim] = sim_lqr(...
            n_samples, n_delay, doff,...
            Ap, Bp, Cp,... % plant
            Aobs, Bobs, Cobs, Kf, Kc,... % observer and regulator
            id_to_bpm, id_to_cm,...
            false);

[y_sim_mode, u_sim_mode] = sim_lqr_mode(...
            n_samples, n_delay, doff,...
            Ap, Bp, Cp,... % plant
            Aoi, Boi, Coi, Kfi, Kci, UR, VR,... % observer and regulator
            id_to_bpm, id_to_cm,...
            false);

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

figure;
subplot(2,2,1); plot(doff(id_to_bpm,:)'); title('Disturbance LQR');
subplot(2,2,2); plot(y_sim(:,id_to_bpm)); title('Output LQR');
subplot(2,2,3); plot(u_sim*scale_u); title('Input LQR');

figure;
subplot(2,2,1); plot(doff(id_to_bpm,:)'); title('Disturbance LQR MODE');
subplot(2,2,2); plot(y_sim_mode(:,id_to_bpm)); title('Output LQR MODE');
subplot(2,2,3); plot(u_sim_mode*scale_u); title('Input LQR MODE');
    
if sim_IMC
    figure;
    subplot(2,2,1); plot(doff(id_to_bpm,:)'); title('Disturbance IMC');
    subplot(2,2,2); plot(y_sim_imc(:,id_to_bpm)); title('Output IMC');
    subplot(2,2,3); plot(u_sim_imc*scale_u); title('Input IMC');
end

% Sensitivity
inds_ss = floor(n_samples/2):n_samples;
[Sfiltered, Sorig, w_Hz] = estim_sensitivity(y_sim(inds_ss,id_to_bpm)', doff(id_to_bpm,inds_ss), length(inds_ss), 1/Ts);
[Sfilteredm, Sorigm, w_Hzm] = estim_sensitivity(y_sim_mode(inds_ss,id_to_bpm)', doff(id_to_bpm,inds_ss), length(inds_ss), 1/Ts);
fname = sprintf('mpc_data_02092022_nd%d.mat',n_delay);
theo_data = load(strrep(fname,'.mat','_analysis.mat'));
if strcmp(pick_direction, 'vertical')
    Stheo = theo_data.Slqr_y;
else
    Stheo = theo_data.Slqr_x;
end
Smin_theo = 20*log10(abs(min(Stheo,[],2)));
Smax_theo = 20*log10(abs(max(Stheo,[],2)));

figure;
subplot(1,2,1);
semilogx(w_Hz,20*log10(abs(Sorig(:,1))),'b'); hold on;
semilogx(theo_data.w_Hz,Smin_theo,'k',theo_data.w_Hz,Smax_theo,'k'); hold on;
semilogx(w_Hz, Sfiltered(:,1),'r');
grid on; title('LQR')

subplot(1,2,2);
semilogx(w_Hz,20*log10(abs(Sorigm(:,1))),'b'); hold on;
semilogx(theo_data.w_Hz,Smin_theo,'k',theo_data.w_Hz,Smax_theo,'k'); hold on;
semilogx(w_Hz, Sfilteredm(:,1),'r');
grid on; title('LQR MODE')

%% IBM
bpm_n = 0;
nfft = 2000;
F_S = 10072;

[ayoff_avg, ~, cyoff_avg] = get_average_psd(doff(id_to_bpm,inds_ss)', length(inds_ss), F_S);
[aysim_avg, ~, cysim_avg] = get_average_psd(y_sim(inds_ss,id_to_bpm), length(inds_ss), F_S);
[aysim_mode_avg, ~, cysim_mode_avg] = get_average_psd(y_sim_mode(inds_ss,id_to_bpm), length(inds_ss), F_S);

lw = 1;
font = 12;

figure1 = figure;
clf;
set(figure1,'Position',[27   300   800   600]);
axes1 = axes('Parent',figure1);
hold(axes1 ,'on');box(axes1 ,'on')
semilogx(ayoff_avg, cyoff_avg,'Parent',axes1,'DisplayName','FOFB OFF','Marker','none',...
    'LineWidth',lw,...
    'LineStyle','-',...
    'Color','r');
semilogx(aysim_avg, cysim_avg,'Parent',axes1,'DisplayName','LQR','Marker','none',...
    'LineWidth',lw,...
    'LineStyle','-',...
    'Color','b');
semilogx(aysim_mode_avg, cysim_mode_avg,'Parent',axes1,'DisplayName','LQR MODE','Marker','none',...
    'LineWidth',lw,...
    'LineStyle','-',...
    'Color','m');
title(sprintf('%s Frequency Response',pick_direction),'FontSize',font,'FontWeight','normal');
ylabel('IBM (\mum)','FontSize',font);
xlabel('Frequency (Hz)','FontSize',font);
grid(axes1,'on');
set(axes1,'XMinorTick','on','Xscale','log');
legend1 = legend(axes1,'show');
set(legend1,'Location','Northwest','FontSize',font-2)
