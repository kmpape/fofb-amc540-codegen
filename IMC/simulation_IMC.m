clc
%close all
clear all
addpath('..')

%% Options
do_codegen = 1; % print files
pick_dir = 1;
dirs = {'X','Y'};
fname_RM = '../ORMS/GoldenBPMResp_I04.mat';
fname_RM_plant = '../ORMS/GoldenBPMResp_I04.mat';
%fname_RM = '../ORMS/GoldenBPMResp_I04.mat';
%fname_RM_plant = '../ORMS/GoldenBPMResp_I04.mat';
codegen_mc = 1; % switch between multicore (mc) or single-core code gen

sr_config_choice = 5;

FULL_CONFIG = 0; % 171 x 172
V3_CONFIG   = 3; % 96 x 96
V4_CONFIG   = 4; % 03.08.2022 SR config
V5_CONFIG   = 5; % 13.09.2022 SR config

%% Load disturbance data
fname_FA = sprintf('../../DATA/26072022_224150_data_%s.mat',dirs{pick_dir}); % FOFB ON
fname_FA = sprintf('../../DATA/03082022_033600_data_%s.mat',dirs{pick_dir}); % FOFB OFF
fname_FA = sprintf('../../DATA/28092022_002500_data_%s.mat',dirs{pick_dir}); % FOFB OFF
ind_start = 1;
ind_end = 100000;

data_FA = load(fname_FA);
if ind_end<0; ind_end = size(data_FA.y, 2); end
inds = ind_start:ind_end;
n_samples = length(inds);
doff = data_FA.y(:,inds);
ufofb = data_FA.u(:,inds);
clear data_FA

%% Configure Diamond-I Storage Ring
% load ORM for controller
load(fname_RM);
RMorigx = Rmat(1).Data(:,:);%  * 1e6; % RM(1) eq to RM(1,1)
[ny_x,nu_x] = size(RMorigx);
RMorigy = Rmat(4).Data(:,:);%  * 1e6; % RM(4) eq to RM(2,2)
[ny_y,nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
[TOT_BPM,TOT_CM] = size(RMorigx);

if (sr_config_choice == V3_CONFIG) % use same storage ring config as GSVD-IMC
    [id_to_bpm_tmp, id_to_cm_tmp, ~] = diamond_I_configuration_v3(RMorigx,RMorigy);
    id_to_bpm_x = id_to_bpm_tmp;
    id_to_bpm_y = id_to_bpm_tmp;
    id_to_cm_x = id_to_cm_tmp;
    id_to_cm_y = id_to_cm_tmp;
    n_cores = 4;
elseif (sr_config_choice == FULL_CONFIG)
    id_to_bpm_tmp = 1:1:173;
    bad_bpm = [76,79];
    id_to_bpm_tmp(bad_bpm) = [];
    id_to_cm_tmp = 1:1:172;
    id_to_bpm_x = id_to_bpm_tmp;
    id_to_bpm_y = id_to_bpm_tmp;
    id_to_cm_x = id_to_cm_tmp;
    id_to_cm_y = id_to_cm_tmp;
    n_cores = 6;
elseif (sr_config_choice == V4_CONFIG) % use same storage ring config as GSVD-IMC
    [id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v4(RMorigx,RMorigy);
    n_cores = 6;
elseif (sr_config_choice == V5_CONFIG) % use same storage ring config as MPC
    square_config = true;
    [id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v5(RMorigx,RMorigy,square_config);
    n_cores = 6;
else
    assert(0);
end

RMx = RMorigx(id_to_bpm_x, id_to_cm_x);
RMy = RMorigy(id_to_bpm_y, id_to_cm_y);

% load ORM for plant
load(fname_RM_plant);
RMorigx_plant = Rmat(1).Data(:,:);%  * 1e6; % RM(1) eq to RM(1,1)
RMorigy_plant = Rmat(4).Data(:,:);%  * 1e6; % RM(4) eq to RM(2,2)

%% Actuators
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time[s]
n_delay = 9; % number of delay time steps [-]
z = tf('z');
s = tf('s');
aIx = 2*pi*500;
aIy = 2*pi*700;
tf_DIx = aIx/(s+aIx);
tf_DIy = aIy/(s+aIy);
gI_mp_zx = c2d(tf_DIx, Ts, 'zoh');
gI_mp_zy = c2d(tf_DIy, Ts, 'zoh');

% add oddities
aslowx = 2*pi*10;
afastx = 2*pi*2000;
gslow_mp_zx = c2d(aslowx/(s+aslowx), Ts, 'zoh');
gfast_mp_zx = c2d(afastx/(s+afastx), Ts, 'zoh');

aslowy = 2*pi*20;
afasty = 2*pi*2000;
gslow_mp_zy = c2d(aslowy/(s+aslowy), Ts, 'zoh');
gfast_mp_zy = c2d(afasty/(s+afasty), Ts, 'zoh');

% input is negated on PMCs, so multiply with -1 here
minus_one = -1;

%% IMC
[U, S, V] = svd(RMx, 'econ');
MU = 1.0*eye(size(RMx,1));
E = S / (S.^2+MU);
Kx = V*E*U';

[U, S, V] = svd(RMy, 'econ');
MU = 1.0*eye(size(RMy,1));
E = S / (S.^2+MU);
Ky = V*E*U';

z = tf('z', Ts);
s = tf('s');
bw = 2*pi*200;
abw = exp(-bw*Ts);
T_mp_z = (1-abw) / (1-z^(-1)*abw) * z^(-1);
q_zx = T_mp_z / gI_mp_zx;
q_zy = T_mp_z / gI_mp_zy;

c_zx = q_zx / (1 - T_mp_z*z^(-n_delay));
c_zy = q_zy / (1 - T_mp_z*z^(-n_delay));

%% Simulation
network_scaling = 1e6;
if pick_dir == 1
    K = Kx / network_scaling;
    cz = c_zx;
    RMorig = RMorigx;
    RMorig_plant = RMorigx_plant;
    gI_mp_z = gI_mp_zx;
    id_to_bpm = id_to_bpm_x;
    id_to_cm = id_to_cm_x;
    gslow_mp_z = gslow_mp_zx;
    gfast_mp_z = gfast_mp_zx;
else
    K = Ky / network_scaling;
    cz = c_zy;
    RMorig = RMorigy;
    RMorig_plant = RMorigy_plant;
    gI_mp_z = gI_mp_zy;
    id_to_bpm = id_to_bpm_y;
    id_to_cm = id_to_cm_y;
    gslow_mp_z = gslow_mp_zy;
    gfast_mp_z = gfast_mp_zy;
end

if true
    [UR,SR,VR] = svd(RMorig_plant);
    bpmdist = find(abs(UR(24,:)) == max(abs(UR(24,:))));
    UR(:,bpmdist) = UR(:,bpmdist)+randn(size(UR(:,bpmdist)))*1e-4;
    UR(:,bpmdist) = UR(:,bpmdist)/norm(UR(:,bpmdist));
    RMtilde = UR*SR*VR';
    [A, B, C ,D] = ssdata([RMorig_plant(:,1:7).*gI_mp_z,...
                           RMorig_plant(:,8:10).*gfast_mp_z,...
                           RMorig_plant(:,11:12).*gslow_mp_z,...
                           RMorig_plant(:,13:172).*gI_mp_z]);
else
    [A, B, C ,D] = ssdata(RMorig .* gI_mp_z);
end

[Ac, Bc, Cc, Dc] = ssdata(-K.*cz);
open_loop = false;
run_ntimes = 1;
[y_sim, u_sim] = sim_standard_imc(...
                    n_samples, n_delay, id_to_bpm, id_to_cm, doff,...
                    Ac, Bc, Cc, Dc,... % CONTROLLER STATE-SPACE
                    A, B, C, D, network_scaling,... % PLANT STATE-SPACE
                    open_loop, 1);
yss = mean(y_sim(end-20:end,:),1);    
t = (0:1:n_samples-1)*Ts;
bad_bpm = [76,79];
all_bpm = 1:1:TOT_BPM; not_ctr = all_bpm; not_ctr([id_to_bpm,bad_bpm])=[];
fig = figure;
subplot(1,4,1); 
if ~isempty(not_ctr); plot(t,doff(not_ctr,:)','--','color',[0.9 0.9 0.9]); hold on; end
plot(t,doff(id_to_bpm,:)'); title(sprintf('Disturbance %s',dirs{pick_dir}));
ylim([min(min(doff(id_to_bpm,:))), max(max(doff(id_to_bpm,:)))]);
xlabel('Time [s]'); ylabel('Position [mum]');

subplot(1,4,2);
if ~isempty(not_ctr); plot(t,y_sim(:,not_ctr),'--','color',[0.9 0.9 0.9]); hold on; end
plot(t,y_sim(:,id_to_bpm)); title(sprintf('Output %s',dirs{pick_dir}))
ylim([min(min(y_sim(:,id_to_bpm))), max(max(y_sim(:,id_to_bpm)))]);
xlabel('Time [s]'); ylabel('Position [mum]');

subplot(1,4,3); plot(not_ctr,yss(not_ctr),'--','color',[0.9 0.9 0.9]); hold on; 
plot(id_to_bpm,yss(id_to_bpm)); title(sprintf('Steady-State %s',dirs{pick_dir}))
ylim([min(min(yss(id_to_bpm))), max(max(yss(id_to_bpm)))]);
xlabel('BPM [-]'); xlim([1 TOT_BPM]); ylabel('Position [mum]');

subplot(1,4,4); plot(t,u_sim); title(sprintf('Setpoint %s',dirs{pick_dir}))
xlabel('Time [s]'); ylabel('Setpoint [A]');

%% PSDs
%cmnt = 'DIAD model I4 plant';
%save('restmp.mat','doff','y_sim','cmnt')
ind_start_fft = 50000;
addpath('../UTILS')
nfft = 10000;
[fs, psd_off, cs_off] = get_all_psd(doff(id_to_bpm,ind_start_fft:end)', nfft, Fs);
[fs, psd_on,  cs_on] = get_all_psd(y_sim(ind_start_fft:end,id_to_bpm), nfft, Fs);

colors = lines(10);
figure; c=1;
semilogx(fs,20*log10(abs(psd_off(40,:))),'Color', colors(c,:),'Linestyle','--'); hold on; c=c+1;
semilogx(fs,20*log10(abs(psd_on(40,:))),'Color', colors(c,:),'Linestyle','-'); hold on; c=c+1;
semilogx(fs,20*log10(abs(psd_off(24,:))),'Color', colors(c,:),'Linestyle','--'); hold on; c=c+1;
semilogx(fs,20*log10(abs(psd_on(24,:))),'Color', colors(c,:),'Linestyle','-'); hold on; c=c+1;
legend('BPM41 OFF','BPM41 ON','BPM24 OFF','BPM24 ON')

%%
colors = lines(10);
for i = [5,10,14]
    figure; c=1;
    semilogx(fs,20*log10(abs(psd_off(i,:))),'Color', colors(c,:),'Linestyle','--'); hold on; c=c+1;
    semilogx(fs,20*log10(abs(psd_on(i,:))),'Color', colors(c,:),'Linestyle','-'); hold on; c=c+1;
    title(sprintf('BPM %d',i));
end

