clc
close all
clear all
addpath('..')

%% Options
do_codegen = 1; % print files
pick_dir = 2;
dirs = {'X','Y'};
fname_RM = '../ORMS/GoldenBPMResp_DIAD.mat';
folder_out = 'out/'; % output folder for codegen
gen_ctr_test_data = 1;
codegen_mc = 1; % switch between multicore (mc) or single-core code gen
gen_filter_unit_test_data = 0;
gen_matvec_mpy_unit_test_data = 0; 
use_PYTHON_parameter = 0;

%% Load disturbance data
fname_FA = sprintf('../../DATA/26072022_224150_data_%s.mat',dirs{pick_dir}); % FOFB OFF
ind_start = 1; 
ind_end = 20000;

data_FA = load(fname_FA);
if ind_end<0; ind_end = size(data_FA.y, 2); end
inds = ind_start:ind_end;
n_samples = length(inds);
doff = data_FA.y(:,inds);
ufofb = data_FA.u(:,inds);
clear data_FA

%% Configure Diamond-I Storage Ring
load(fname_RM);
RMorigx = Rmat(1).Data(:,:);%  * 1e6; % RM(1) eq to RM(1,1)
[ny_x,nu_x] = size(RMorigx);
RMorigy = Rmat(4).Data(:,:);%  * 1e6; % RM(4) eq to RM(2,2)
[ny_y,nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
[TOT_BPM,TOT_CM] = size(RMorigx);
if (0) % use same storage ring config as GSVD-IMC
    [id_to_bpm, slow_to_id, ~] = diamond_I_configuration_v3(RMorigx,RMorigy);
    n_cores = 4;
else % full storage ring
    id_to_bpm = 1:1:173;
    bad_bpm = [76,79];
    id_to_bpm(bad_bpm) = [];
    slow_to_id = 1:1:172;
    n_cores = 6;
end

RMx = RMorigx(id_to_bpm, slow_to_id);
RMy = RMorigy(id_to_bpm, slow_to_id);

ny = length(id_to_bpm);
nu = length(slow_to_id);

%% Actuators
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time[s]
n_delay = 8; % number of delay time steps [-]
z = tf('z');
s = tf('s');
aIx = 2*pi*500;
aIy = 2*pi*700;
tf_DIx = aIx/(s+aIx);
tf_DIy = aIy/(s+aIy);
gI_mp_zx = c2d(tf_DIx, Ts, 'zoh');
gI_mp_zy = c2d(tf_DIy, Ts, 'zoh');

% input is negated on PMCs, so multiply with -1 here
minus_one = -1;

%% IMC
[U, S, V] = svd(RMx, 'econ');
MU = 1.0*eye(ny);
E = S / (S.^2+MU);
Kx = V*E*U';

[U, S, V] = svd(RMy, 'econ');
MU = 1.0*eye(ny);
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
    gI_mp_z = gI_mp_zx;
else
    K = Ky / network_scaling;
    cz = c_zy;
    RMorig = RMorigy;
    gI_mp_z = gI_mp_zy;
end

[A, B, C ,D] = ssdata(RMorig .* gI_mp_z);
[Ac, Bc, Cc, Dc] = ssdata(-K.*cz);
[y_sim, u_sim] = sim_standard_imc(...
                    n_samples, n_delay, id_to_bpm, slow_to_id, doff,...
                    Ac, Bc, Cc, Dc,... % CONTROLLER STATE-SPACE
                    A, B, C, D, network_scaling); % PLANT STATE-SPACE

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



