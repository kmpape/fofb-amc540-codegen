clc
close all
clear all
addpath('..')

%% Options
do_codegen = 1; % print files
fname_RM = '../ORMS/ORM_11.9.2022/GoldenBPMResp_I04.mat';
folder_out = 'out/'; % output folder for codegen
codegen_mc = 1; % switch between multicore (mc) or single-core code gen
gen_filter_unit_test_data = 0;
gen_matvec_mpy_unit_test_data = 0;
gen_ctr_test_data = 1;
use_PYTHON_parameter = 0;

%% Configure Diamond-I Storage Ring
load(fname_RM);
RMorigx = Rmat(1).Data(:,:);%  * 1e6; % RM(1) eq to RM(1,1)
[ny_x, nu_x] = size(RMorigx);
RMorigy = Rmat(4).Data(:,:);%  * 1e6; % RM(4) eq to RM(2,2)
[ny_y, nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
[TOT_BPM, TOT_CM] = size(RMorigx);
[id_to_bpm, slow_to_id, fast_to_id] = diamond_I_configuration_v3(RMorigx,RMorigy);
Rsx = RMorigx(id_to_bpm, slow_to_id);
Rsy = RMorigy(id_to_bpm, slow_to_id);
Rfx = RMorigx(id_to_bpm, fast_to_id);
Rfy = RMorigy(id_to_bpm, fast_to_id);

ny = length(id_to_bpm);
ns = length(slow_to_id);
nf = length(fast_to_id);

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

% input is negated on PMCs, so multiply with -1 here
minus_one = -1;

%% Mid-Ranging IMC
bw_allx = 1/(n_delay*Ts);%176*2*pi; % overall desired bandwidth [rad/s]
bw_sx =  10*2*pi; % slow actuators desired bandwidth [rad/s]
T_tiso_mpx = bw_allx/(s+bw_allx);
T_siso_mpx = bw_sx/(s+bw_sx);

bw_ally = 1/(n_delay*Ts);%176*2*pi; % overall desired bandwidth [rad/s]
bw_sy = 10*2*pi; % slow actuators desired bandwidth [rad/s]
T_tiso_mpy = bw_ally/(s+bw_ally);
T_siso_mpy = bw_sy/(s+bw_sy);

% new: epsilon parameter as in the standard FOFB
epsilon = 1e-6;
qs_zx = (1-epsilon) * c2d(T_siso_mpx / tf_DIx, Ts, 'zoh');
qf_zx = (1-epsilon) * c2d((T_tiso_mpx - T_siso_mpx) / tf_DIx, Ts, 'zoh');
qs_zy = (1-epsilon) * c2d(T_siso_mpy / tf_DIy, Ts, 'zoh');
qf_zy = (1-epsilon) * c2d((T_tiso_mpy - T_siso_mpy) / tf_DIy, Ts, 'zoh');

% GSVD
mu = 1;

[Us,Uf,X,C,S] = gsvd(Rsx', Rfx');
Ss = C';
Sf = S';
Ff = X*pinv(X*blkdiag(eye(rank(Rfx)), zeros(ns-rank(Rfx))));
G = X*((X'*X+mu*eye(ny))\X');
Ksx = -(Us*inv(Ss)/X)*G;
Kfx = -((Uf*pinv(Sf)/X)*Ff)*G;
Psx = -G\Rsx;
Pfx = -G\Rfx;

[Us,Uf,X,C,S] = gsvd(Rsy', Rfy');
Ss = C';
Sf = S';
Ff = X*pinv(X*blkdiag(eye(rank(Rfy)), zeros(ns-rank(Rfy))));
G = X*((X'*X+mu*eye(ny))\X');
Ksy = -(Us*inv(Ss)/X)*G;
Kfy = -((Uf*pinv(Sf)/X)*Ff)*G;
Psy = -G\Rsy;
Pfy = -G\Rfy;

w_Hz = logspace(-2,3,200);
get_GSVD_IMC_sensitivity_v2(RMorigx, RMorigy, bw_sx, bw_allx, bw_sy, bw_ally, n_delay, w_Hz, true, n_delay, mu);

%% Code generation
n_cores = 4;
if codegen_mc == true
    ns_per_core = 32;
    ny_per_core = 32;
    nf_per_core = 16;
else
    ns_per_core = n_cores*32;
    ny_per_core = n_cores*32;
    nf_per_core = n_cores*16;
end
ns_pad = n_cores*32-ns;
nf_pad = n_cores*16-nf;
ny_pad = ns_pad;

assert(ns_pad >= 0);
assert(nf_pad >= 0);

dirs = {'x','y'};
filter_len = ns_per_core;
network_scaling = 1e6;
if do_codegen == true
    [num_qsx, den_qsx] = tfdata(qs_zx);
    [num_qsy, den_qsy] = tfdata(qs_zy);
    generate_filter_from_TF_XY(folder_out, filter_len,...
        num_qsx{:}, den_qsx{:}, 'qs', 'double', gen_filter_unit_test_data,...
        '.gsvd_qs', 64, 4.99, num_qsy{:}, den_qsy{:});
    
    [num_qfx, den_qfx] = tfdata(qf_zx);
    [num_qfy, den_qfy] = tfdata(qf_zy);
    generate_filter_from_TF_XY(folder_out, filter_len,...
        num_qfx{:}, den_qfx{:}, 'qf', 'double', gen_filter_unit_test_data,...
        '.gsvd_qf', 64, 4.99, num_qfy{:}, den_qfy{:});
    
    % NOTE: for D-I slow and fast models are identical
    [num_gsx, den_gsx] = tfdata(gI_mp_zx .* z^(-n_delay));
    [num_gsy, den_gsy] = tfdata(gI_mp_zy .* z^(-n_delay));
    generate_filter_from_TF_XY(folder_out, filter_len,...
        num_gsx{:}, den_gsx{:}, 'gs', 'double', gen_filter_unit_test_data,...
        '.gsvd_gs', 64, 1e20, num_gsy{:}, den_gsy{:});
    
    [num_gfx, den_gfx] = tfdata(gI_mp_zx .* z^(-n_delay));
    [num_gfy, den_gfy] = tfdata(gI_mp_zy .* z^(-n_delay));
    generate_filter_from_TF_XY(folder_out, filter_len,...
        num_gfx{:}, den_gfx{:}, 'gf', 'double', gen_filter_unit_test_data,...
        '.gsvd_gf', 64, 1e20, num_gfy{:}, den_gfy{:});
    
    for i = 1 : 2
        
        % print gain matrix
        fid = fopen([folder_out,'GSVD_gain_',dirs{i},'.h'], 'w');
        fprintf(fid, '#ifndef GSVD_GAIN_H\n#define GSVD_GAIN_H\n');
        fprintf(fid, '#define GSVD_GAIN_NY (%d)\n', TOT_BPM);
        fprintf(fid, '#define GSVD_GAIN_NU (%d)\n', TOT_CM);
        fprintf(fid, '#define GSVD_GAIN_NS (%d)\n', ns);
        fprintf(fid, '#define GSVD_GAIN_NF (%d)\n', nf);
        fprintf(fid, '#define GSVD_GAIN_NS_PAD (%d)\n', ns+ns_pad);
        fprintf(fid, '#define GSVD_GAIN_NF_PAD (%d)\n', nf+nf_pad);

        if i == 1
            Ks = minus_one * blkdiag(Ksx, zeros(ns_pad,ny_pad)) / network_scaling;
            Kf = minus_one * blkdiag(Kfx, zeros(nf_pad,ny_pad)) / network_scaling;
            Ps = minus_one * blkdiag(Psx, zeros(ny_pad,ns_pad)) * network_scaling;
            Pf = minus_one * blkdiag(Pfx, zeros(ny_pad,nf_pad)) * network_scaling;
        else
            Ks = minus_one * blkdiag(Ksy, zeros(ns_pad,ny_pad)) / network_scaling;
            Kf = minus_one * blkdiag(Kfy, zeros(nf_pad,ny_pad)) / network_scaling;
            Ps = minus_one * blkdiag(Psy, zeros(ny_pad,ns_pad)) * network_scaling;
            Pf = minus_one * blkdiag(Pfy, zeros(ny_pad,nf_pad)) * network_scaling;
        end
        
        % NOTE:
        % (1) PMCs multiplies the signal passed to the magnets by -1 (see
        % minus_one above)
        % (2) Input to the magnets is in [A]
        % (3) BPM signals are passed in [nm]
        % Here ORM in [nm/A]
        print_dense_C_matrix(fid, Ks, 'float', 'GSVD_Ks_gain', true,'.gsvd_init', 2);
        print_dense_C_matrix(fid, Kf, 'float', 'GSVD_Kf_gain', true,'.gsvd_init', 2);
        print_dense_C_matrix(fid, Ps, 'float', 'GSVD_Ps_gain', true,'.gsvd_init', 2);
        print_dense_C_matrix(fid, Pf, 'float', 'GSVD_Pf_gain', true,'.gsvd_init', 2);
        fprintf(fid, '#endif\n');
        fclose(fid);
        
        fid = fopen([folder_out,'GSVD_storage_ring_config_',dirs{i},'.h'], 'w');
        fprintf(fid, '#ifndef GSVD_STORAGE_RING_CONFIG_H\n#define GSVD_STORAGE_RING_CONFIG_H\n');
        print_dense_C_matrix(fid, id_to_bpm-1, 'int', 'GSVD_ID_TO_BPM', true,'.gsvd_shared', 2);
        print_dense_C_matrix(fid, slow_to_id-1, 'int', 'GSVD_SLOW_TO_BPM', true,'.gsvd_shared', 2);
        print_dense_C_matrix(fid, fast_to_id-1, 'int', 'GSVD_FAST_TO_BPM', true,'.gsvd_shared', 2);
        fprintf(fid, '#endif\n');
        fclose(fid);
        
        if gen_matvec_mpy_unit_test_data == 1
            fid = fopen([folder_out,'GSVD_UNIT_TEST_',dirs{i},'.h'], 'w');
            fprintf(fid, '#ifndef GSVD_UNIT_TEST_H\n#define GSVD_UNIT_TEST_H\n');
            test_in = randn(size(Ks,2),1); test_out = Ks * test_in;
            print_dense_C_matrix(fid, test_in, 'double', 'GSVD_test_in_Ks', true,'.gsvd_unit_test', 2);
            print_dense_C_matrix(fid, test_out, 'double', 'GSVD_test_out_Ks', true,'.gsvd_unit_test', 2);
            
            test_in = randn(size(Kf,2),1); test_out = Kf * test_in;
            print_dense_C_matrix(fid, test_in, 'double', 'GSVD_test_in_Kf', true,'.gsvd_unit_test', 2);
            print_dense_C_matrix(fid, test_out, 'double', 'GSVD_test_out_Kf', true,'.gsvd_unit_test', 2);
            
            test_in = randn(size(Ps,2),1); test_out = Ps * test_in;
            print_dense_C_matrix(fid, test_in, 'double', 'GSVD_test_in_Ps', true,'.gsvd_unit_test', 2);
            print_dense_C_matrix(fid, test_out, 'double', 'GSVD_test_out_Ps', true,'.gsvd_unit_test', 2);
            
            test_in = randn(size(Pf,2),1); test_out = Pf * test_in;
            print_dense_C_matrix(fid, test_in, 'double', 'GSVD_test_in_Pf', true,'.gsvd_unit_test', 2);
            print_dense_C_matrix(fid, test_out, 'double', 'GSVD_test_out_Pf', true,'.gsvd_unit_test', 2);
            
            fprintf(fid, '#endif\n');
        end
        
        if gen_ctr_test_data == 1
            n_samples = 1000; 
            %doff = randn(TOT_BPM,1).*ones(1,n_samples);
            doff = ones(TOT_BPM,1).*ones(1,n_samples);
            if i==1
                tf_DI = tf_DIx; 
                RMorig = RMorigx;
                T_tiso_mp = T_tiso_mpx;
                T_siso_mp = T_siso_mpx;
            else
                tf_DI = tf_DIy; 
                RMorig = RMorigy;
                T_tiso_mp = T_tiso_mpy;
                T_siso_mp = T_siso_mpy;
            end
            gI_mp_z = c2d(tf_DI, Ts, 'zoh');
            [A, B, C ,D] = ssdata(RMorig .* gI_mp_z);
            [Acs,Bcs,Ccs,Dcs] = ssdata(eye(ns) .* c2d(T_siso_mp / tf_DI, Ts, 'zoh'));
            [Acf,Bcf,Ccf,Dcf] = ssdata(eye(nf) .* c2d((T_tiso_mp - T_siso_mp) / tf_DI, Ts, 'zoh'));
            [Ams,Bms,Cms,Dms] = ssdata(eye(ny) .* gI_mp_z .* z^(-n_delay));
            [Amf,Bmf,Cmf,Dmf] = ssdata(eye(ny) .* gI_mp_z .* z^(-n_delay));

            [y_sim, u_sim, us_sim, uf_sim] = sim_mid_range_gsvd(...
                n_delay, doff,...
                A, B, C, D,...               % Plant
                -Ks(1:end-ns_pad,1:end-ny_pad), Acs, Bcs, Ccs, Dcs,...       % Slow controller
                -Kf(1:end-nf_pad,1:end-ny_pad), Acf, Bcf, Ccf, Dcf,...       % Fast controller
                -Ps(1:end-ny_pad,1:end-ns_pad), Ams, Bms, Cms, Dms,...       % Slow plant model
                -Pf(1:end-ny_pad,1:end-nf_pad), Amf, Bmf, Cmf, Dmf,...       % Fast plant model
                id_to_bpm, slow_to_id, fast_to_id, false, network_scaling);
            nt = 100;
            ydata = round(y_sim(1:1+nt,:)*1000,0);
            udata = round(u_sim(2+n_delay:2+n_delay+nt,:)*network_scaling,0);
            fid = fopen([folder_out,'GSVD_test_data_',dirs{i},'.h'], 'w');
            fprintf(fid, sprintf('#ifndef GSVD_TEST_DATA_%s_H\n#define GSVD_TEST_DATA_%s_H\n',upper(dirs{i}),upper(dirs{i})));
            fprintf(fid, '#define GSVD_NTEST (%d)\n', nt);
            print_dense_C_matrix(fid, ydata, 'int', 'GSVD_TEST_IN', true,'.gsvd_unit_test', 2);
            print_dense_C_matrix(fid, udata, 'int', 'GSVD_TEST_OUT', true,'.gsvd_unit_test', 2);
            fprintf(fid, '#endif\n');
            fclose(fid);
            
            yss = mean(y_sim(end-20:end,:),1);    
            t = (0:1:n_samples-1)*Ts;
            bad_bpm = [76,79];
            all_bpm = 1:1:TOT_BPM; not_ctr = all_bpm; not_ctr([id_to_bpm,bad_bpm])=[];
            fig = figure;
            subplot(1,5,1); plot(t,doff(not_ctr,:)','--','color',[0.9 0.9 0.9]); hold on; 
            plot(t,doff(id_to_bpm,:)'); title(sprintf('Disturbance %s',dirs{i}));
            ylim([0.9*min(min(doff(id_to_bpm,:))), 1.1*max(max(doff(id_to_bpm,:)))]);
            xlabel('Time [s]'); ylabel('Position [mum]');

            subplot(1,5,2); plot(t,y_sim(:,not_ctr),'--','color',[0.9 0.9 0.9]); hold on; 
            plot(t,y_sim(:,id_to_bpm)); title(sprintf('Output %s',dirs{i}))
            ylim([min(min(y_sim(:,id_to_bpm))), max(max(y_sim(:,id_to_bpm)))]);
            xlabel('Time [s]'); ylabel('Position [mum]');

            subplot(1,5,3); plot(not_ctr,yss(not_ctr),'--','color',[0.9 0.9 0.9]); hold on; 
            plot(id_to_bpm,yss(id_to_bpm)); title(sprintf('Steady-State %s',dirs{i}))
            ylim([min(min(yss(id_to_bpm))), max(max(yss(id_to_bpm)))]);
            xlabel('BPM [-]'); xlim([1 TOT_BPM]); ylabel('Position [mum]');

            subplot(1,5,4); plot(t,us_sim); title(sprintf('Slow Setpoint %s',dirs{i}))
            xlabel('Time [s]'); ylabel('Setpoint [A]');
            subplot(1,5,5); plot(t,uf_sim); title(sprintf('Fast Setpoint %s',dirs{i}))
            xlabel('Time [s]'); ylabel('Setpoint [A]');
        end
    end
end



