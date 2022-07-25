clc
%close all
clear all
addpath('..')

%% Options
do_codegen = 1; % print files
fname_RM = '../ORMS/GoldenBPMResp_DIAD.mat';
folder_out = 'out/'; % output folder for codegen
gen_ctr_test_data = 1;
codegen_mc = 1; % switch between multicore (mc) or single-core code gen
gen_filter_unit_test_data = 0;
gen_matvec_mpy_unit_test_data = 0; 
use_PYTHON_parameter = 0;

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
else % full storage ring
    id_to_bpm = 1:1:173;
    bad_bpm = [76,79];
    id_to_bpm(bad_bpm) = [];
    slow_to_id = 1:1:172;
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

%% Code generation
n_cores = 6;
if codegen_mc == true
    nu_per_core = 32;
    ny_per_core = 32;
else
    nu_per_core = n_cores*32;
    ny_per_core = n_cores*32;
end
nu_pad = n_cores*32-nu;
ny_pad = n_cores*32-ny;

dirs = {'x','y'};
network_scaling = 1e6;
if do_codegen == true
    for i = 1 : 2
        [num_czx, den_czx] = tfdata(c_zx);
        [num_czy, den_czy] = tfdata(c_zy);
        % generate filter
        generate_filter_from_TF_XY(folder_out, nu, num_czx{:}, den_czx{:}, 'IMC_DI',...
            'double', gen_filter_unit_test_data, '.imc_DI', 64, 4.99, num_czy{:}, den_czy{:});
        % generate_filter_from_TF_XY(foldername, vec_len, num_f, den_f, filter_name,...
             % filter_base_type, print_test_data, datasection, alignment, integrator_max, num_f2, den_f2)
        % print gain matrix
        if codegen_mc == true
            fid = fopen([folder_out,'IMC_DI_gain_mc_',dirs{i},'.h'], 'w');
        else
            fid = fopen([folder_out,'IMC_DI_gain_',dirs{i},'.h'], 'w');
        end
        fprintf(fid, '#ifndef IMC_DI_GAIN_H\n#define IMC_DI_GAIN_H\n');
        fprintf(fid, '#define IMC_DI_GAIN_NY (%d)\n', ny);
        fprintf(fid, '#define IMC_DI_GAIN_NU (%d)\n', nu);
        
        if i == 1
            K = blkdiag(Kx, zeros(nu_pad,ny_pad)) / network_scaling;
        else
            K = blkdiag(Ky, zeros(nu_pad,ny_pad)) / network_scaling;
        end
        
        % NOTE:
        % (1) PMCs multiplies the signal passed to the magnets by -1
        % (2) Input to the magnets is in [A]
        % (3) BPM signals are passed in [nm]
        % Here ORM in [nm/A]
        print_dense_C_matrix(fid, K, 'float', 'IMC_DI_gain_mat', true, '.imc_DI_init', 2);
        fprintf(fid, '#endif\n');
        fclose(fid);
        
        fid = fopen([folder_out,'IMC_storage_ring_config_',dirs{i},'.h'], 'w');
        fprintf(fid, '#ifndef IMC_STORAGE_RING_CONFIG_H\n#define IMC_STORAGE_RING_CONFIG_H\n');
        print_dense_C_matrix(fid, id_to_bpm-1, 'int', 'IMC_ID_TO_BPM', true, '.imc_shared', 2);
        print_dense_C_matrix(fid, slow_to_id-1, 'int', 'IMC_CM_TO_BPM', true, '.imc_shared', 2);
        fprintf(fid, '#endif\n');
        fclose(fid);
        
        if (gen_ctr_test_data == 1)
            if i == 1
                cz = c_zx;
                RMorig = RMorigx;
                gI_mp_z = gI_mp_zx;
            else
                cz = c_zy;
                RMorig = RMorigy;
                gI_mp_z = gI_mp_zy;
            end
            n_samples = 1000; 
            doff = randn(TOT_BPM,1).*ones(1,n_samples);
            [A, B, C ,D] = ssdata(RMorig .* gI_mp_z);
            [Ac, Bc, Cc, Dc] = ssdata(-K(1:end-nu_pad,1:end-ny_pad).*cz);
            [y_sim, u_sim] = sim_standard_imc(...
                                n_samples, n_delay, id_to_bpm, slow_to_id, doff,...
                                Ac, Bc, Cc, Dc,... % CONTROLLER STATE-SPACE
                                A, B, C, D, network_scaling); % PLANT STATE-SPACE
                    
            nt = 50;
            ydata = round(y_sim(1:1+nt,:)*1000,0);
            udata = round(u_sim(9:9+nt,:)*1e6,0);
            fid = fopen([folder_out,'IMC_test_data_',dirs{i},'.h'], 'w');
            fprintf(fid, sprintf('#ifndef IMC_TEST_DATA_%s_H\n#define IMC_TEST_DATA_%s_H\n',upper(dirs{i}),upper(dirs{i})));
            fprintf(fid, '#define IMC_NTEST (%d)\n', nt);
            print_dense_C_matrix(fid, ydata, 'int', 'IMC_TEST_IN', true,'.imc_unit_test', 2);
            print_dense_C_matrix(fid, udata, 'int', 'IMC_TEST_OUT', true,'.imc_unit_test', 2);
            fprintf(fid, '#endif\n');
            fclose(fid);
            
            yss = mean(y_sim(end-20:end,:),1);    
            t = (0:1:n_samples-1)*Ts;
            bad_bpm = [76,79];
            all_bpm = 1:1:TOT_BPM; not_ctr = all_bpm; not_ctr([id_to_bpm,bad_bpm])=[];
            fig = figure;
            subplot(1,4,1); 
            if ~isempty(not_ctr); plot(t,doff(not_ctr,:)','--','color',[0.9 0.9 0.9]); hold on; end
            plot(t,doff(id_to_bpm,:)'); title(sprintf('Disturbance %s',dirs{i}));
            ylim([min(min(doff(id_to_bpm,:))), max(max(doff(id_to_bpm,:)))]);
            xlabel('Time [s]'); ylabel('Position [mum]');

            subplot(1,4,2);
            if ~isempty(not_ctr); plot(t,y_sim(:,not_ctr),'--','color',[0.9 0.9 0.9]); hold on; end
            plot(t,y_sim(:,id_to_bpm)); title(sprintf('Output %s',dirs{i}))
            ylim([min(min(y_sim(:,id_to_bpm))), max(max(y_sim(:,id_to_bpm)))]);
            xlabel('Time [s]'); ylabel('Position [mum]');

            subplot(1,4,3); plot(not_ctr,yss(not_ctr),'--','color',[0.9 0.9 0.9]); hold on; 
            plot(id_to_bpm,yss(id_to_bpm)); title(sprintf('Steady-State %s',dirs{i}))
            ylim([min(min(yss(id_to_bpm))), max(max(yss(id_to_bpm)))]);
            xlabel('BPM [-]'); xlim([1 TOT_BPM]); ylabel('Position [mum]');

            subplot(1,4,4); plot(t,u_sim); title(sprintf('Setpoint %s',dirs{i}))
            xlabel('Time [s]'); ylabel('Setpoint [A]');
        end
    end
end



