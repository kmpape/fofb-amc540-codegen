clc
%close all
clear all
addpath('..')

%% Options
do_codegen = 1; % print files
fname_RM = '../ORMS/GoldenBPMResp_DIAD.mat';
folder_out = 'out/'; % output folder for codegen
n_delay = 9;

%% Hardlimits
load('../ORMS/correctors.mat');
hardlimits = corrector_data.MaxAmps(1:172); % in Amperes

% fid = fopen(sprintf('%s/watchdog_data.h',folder_out), 'w');
% fprintf(fid, '#ifndef WATCHDOG_DATA_H_\n#define WATCHDOG_DATA_H_\n\n');
% print_dense_C_matrix(fid, [hardlimits*1e6; zeros(20,1)], 'int',...
%     'HARDLIMITS_SCALED', true, '.local_data', 2);
% fprintf(fid, '#endif /* WATCHDOG_DATA_H_ */\n');
% fclose(fid);

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

% NOTE: storage configuration, gains etc. must match the data used for
% preconditioning the Hessian (done in precond_LMI_v2.m)

%% Observer and Regulator
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time [s]
fname = sprintf('mpc_data_02092022_nd%d.mat',n_delay);
if ~exist(fname,'file')
    print_msg = false;
    [Ao_x, Bo_x, Co_x, Ap_x, Bp_x, Cp_x, Ad_x, Cd_x,...
     Kfd_x, Kfx_x, ~, ~, P_x, Rlqr_x, Qlqr_x,...
     Ao_y, Bo_y, Co_y, Ap_y, Bp_y, Cp_y, Ad_y, Cd_y,...
     Kfd_y, Kfx_y, ~, ~, P_y, Rlqr_y, Qlqr_y] =...
          observer_regulator(RMorigx,RMorigy,id_to_bpm_x,id_to_cm_x,id_to_bpm_y,id_to_cm_y,n_delay,fname,print_msg);
else
    load(fname);
end

%% Hard-coded data
fprintf('\n#define OBS_delay           (%d) // HARDCODED. CAN BE 8 or 9.\n', n_delay)
fprintf('\n#if (XDIR == 1)\n')
fprintf('const obs_float OBS_Ax = %.16f;\n',Ao_x(1,1))
fprintf('const obs_float OBS_Bx = %.16f;\n',Bo_x(1,1))
for j = 1:9
    fprintf('const obs_float OBS_Ax_pow_%d = %.16f;\n', j, Ao_x(1,1)^j)
end
fprintf('#else\n')
fprintf('const obs_float OBS_Ax = %.16f;\n',Ao_y(1,1))
fprintf('const obs_float OBS_Bx = %.16f;\n',Bo_y(1,1))
for j = 2:9
    fprintf('const obs_float OBS_Ax_pow_%d = %.16f;\n', j, Ao_y(1,1)^j)
end
fprintf('#endif\n')

%%
pick_dirs = {'horizontal','vertical'};
dirs = {'X','Y'};
for pick_dir = 1:2
    pick_direction = pick_dirs{pick_dir};
    
    % NOTE: SOFB setpoints are read from the FPGA
    
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
        
        SOFB_setpoints = 0*1000*[0.017018109560012817, -0.36062681674957275, 0.4265265464782715, -1.2430311441421509, 1.84929358959198, -1.3640027046203613, 1.089625597000122, -1.152754545211792, -3.0912649631500244, 1.712203025817871, -0.5711809396743774, -0.10924321413040161, -0.53889000415802, 0.7232564091682434, -0.4079737663269043, -0.6608858108520508, 0.6217878460884094, -2.5586791038513184, 1.3671751022338867, -0.18846715986728668, -0.19866496324539185, 0.3715630769729614, 0.6067423820495605, -0.7344881296157837, 1.216092824935913, 0.41330069303512573, -0.8586440682411194, -0.43454259634017944, -0.1075853556394577, 0.13515017926692963, -0.04713757336139679, -0.3805052638053894, 0.12361033260822296, 0.008357536047697067, -0.3025299906730652, 0.4496183693408966, -0.4757261574268341, -0.17698454856872559, -0.6850622296333313, 0.8487122654914856, -0.2911084294319153, -1.4305939674377441, 1.1703999042510986, -0.7916262149810791, 0.7223260998725891, 0.9766864776611328, -1.0604008436203003, 0.2158363163471222, -0.7098335027694702, 1.3440566062927246, 0.10605031996965408, -0.13034237921237946, -1.3893193006515503, -0.6569163203239441, 2.5632588863372803, -0.5431381464004517, -0.23142468929290771, 0.10426411777734756, -0.3864889442920685, 1.0281198024749756, 0.660631000995636, -2.4857847690582275, 0.531500518321991, -0.32756051421165466, 2.135746955871582, -5, 0.4553358256816864, -0.6173129081726074, 0.23832263052463531, -0.07054101675748825, 0.07058994472026825, -0.6087528467178345, 0.8179794549942017, 0.4336724281311035, -0.7508447766304016, 1.5841120481491089, -1.0532056093215942, -0.25923895835876465, 0.23892439901828766, 0.38093677163124084, -0.3947981595993042, 0.8433530926704407, 0.03773653507232666, -1.1912553310394287, 0.6339174509048462, -0.21830596029758453, 0.3377736210823059, -0.3769873082637787, 0.3557036519050598, -1.2388724088668823, 1.9396175146102905, -1.0390145778656006, 0.5461440086364746, 0.06459860503673553, -1.4476672410964966, 0.1209847703576088, -0.05696037411689758, 0.3684808313846588, -1.1397080421447754, 2.6863250732421875, -0.41714945435523987, -0.16479560732841492, 0.6453186869621277, -0.5271632075309753, 0.7578499913215637, -1.2476438283920288, 0.6588857173919678, 0.5793570876121521, -0.38712066411972046, -1.0227831602096558, 1.1069272756576538, -1.5502712726593018, 0.39399656653404236, 1.0665626525878906, -0.37509265542030334, -0.2444313019514084, 0.2607043385505676, 0.2630053758621216, 0.7709961533546448, -0.1745031476020813, -0.4268176257610321, -0.5004496574401855, 0.25761422514915466, -0.15045979619026184, -0.3884647488594055, 0.36230260133743286, 0.6295022368431091, -1.0143496990203857, 0.6934900879859924, -0.461463063955307, -0.3647838532924652, 0.18088854849338531, 0.2159477174282074, -0.03209175169467926, -0.42642608284950256, 0.4884761571884155, -0.06146499514579773, -0.08903351426124573, -0.014827894046902657, 0.24346567690372467, -0.43520358204841614, 1.1171962022781372, 0.10287937521934509, -0.2896258533000946, -0.27106067538261414, 0.2312719225883484, 1.637657642364502, -2.497133731842041, 2.1709370613098145, -0.5073715448379517, 0.34836217761039734, 0.6982198357582092, -0.96039879322052, -0.08269193023443222, 0.7715268731117249, -2.186379909515381, 0.6306256055831909, 0.21947143971920013, 0.012526738457381725, -0.38719552755355835, 0.5969539284706116, -0.1951894909143448, 0.024609176442027092, 0.12410219758749008, -0.37557801604270935, -1.2077199220657349, 1.1120285987854004, -0.5590432286262512, -0.93660569190979, 1.9882194995880127, 0.15305931866168976, -0.6962593793869019];
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
        
        SOFB_setpoints = 0*1000*[0.01619904302060604, -0.3885023295879364, 0.6992301344871521, -1.7483879327774048, 2.3418004512786865, -1.5496907234191895, 1.1642824411392212, -1.0826029777526855, -3.5085201263427734, 2.117027997970581, -0.7213708162307739, -0.016367638483643532, -0.8390787839889526, 0.9059506058692932, -0.44695019721984863, -0.6361433267593384, 0.5771667957305908, -2.4178178310394287, 1.3213233947753906, -0.41284069418907166, -0.01976143568754196, 0.28217950463294983, 0.617164134979248, -0.7120676636695862, 1.0387179851531982, 0.4763958752155304, -0.5612682104110718, -0.6951331496238708, 0.02955007180571556, 0.11709681898355484, -0.044195499271154404, -0.45408889651298523, 0.38928934931755066, -0.3819846510887146, -0.11148449033498764, 0.3617706894874573, -0.4519191384315491, -0.20117510855197906, -0.6386896371841431, 0.7802280187606812, -0.22534973919391632, -1.4564428329467773, 1.1819654703140259, -0.7721765637397766, 0.6452361941337585, 1.4174635410308838, -1.711761713027954, 0.6994525194168091, -0.8488370776176453, 1.4066176414489746, 0.07518509030342102, -0.070978082716465, -1.5602132081985474, -0.5684106945991516, 2.698157548904419, -0.6943356990814209, -0.13617253303527832, 0.0965278148651123, -0.4131055474281311, 1.202785849571228, 0.6013646721839905, -2.7761335372924805, 0.7537578344345093, -0.4265475571155548, 2.0365219116210938, -5, 0.4749950170516968, -0.6516678929328918, 0.4227065145969391, -0.342603862285614, 0.27613964676856995, -0.6612204313278198, 0.8361185193061829, 0.4531811773777008, -0.8470264673233032, 1.780594825744629, -1.05325448513031, -0.3227340281009674, 0.2727269232273102, 0.3892798125743866, -0.42446842789649963, 0.97515869140625, -0.010264929383993149, -1.4114340543746948, 0.8374415040016174, -0.3478398621082306, 0.3624727129936218, -0.4587855637073517, 0.7216208577156067, -1.558065414428711, 1.8873403072357178, -0.8834448456764221, 0.4664323329925537, 0.13014689087867737, -1.6719253063201904, 0.13990825414657593, -0.062160685658454895, 0.30241018533706665, -1.0267870426177979, 2.59081768989563, -0.38700228929519653, -0.176784947514534, 0.65232914686203, -0.5602529644966125, 0.9080629944801331, -1.3520822525024414, 0.5552564263343811, 0.6984527707099915, -0.4484162926673889, -0.9921460747718811, 1.0305604934692383, -1.0470114946365356, -0.5739063620567322, 2.045651912689209, -0.7531506419181824, -0.05252191796898842, 0.23519375920295715, 0.2581924498081207, 0.7864719033241272, -0.04462754726409912, -0.7147292494773865, -0.3432796895503998, 0.1804860681295395, -0.1353411078453064, -0.3737372159957886, 0.03207216039299965, 1.3494750261306763, -1.823861002922058, 1.0277291536331177, -0.6163777112960815, -0.30168670415878296, 0.0707305371761322, 0.496951699256897, -0.13955195248126984, -0.736243486404419, 0.7516483068466187, -0.19038543105125427, -0.06523245573043823, -0.02522803097963333, 0.26334622502326965, -0.4521408975124359, 1.0672394037246704, 0.1539207249879837, -0.3190993666648865, -0.2580198645591736, 0.21354863047599792, 1.6735446453094482, -2.5382795333862305, 2.207550287246704, -0.5185757875442505, 0.3553980588912964, 0.691810131072998, -0.9267514944076538, -0.46967199444770813, 1.5809837579727173, -3.0611188411712646, 0.9965386390686035, 0.07120700925588608, 0.05444963276386261, -0.4420344829559326, 0.759113609790802, -0.33723288774490356, -0.0011737275635823607, 0.19979232549667358, -0.414171427488327, -1.1921398639678955, 1.0765595436096191, -0.31856679916381836, -1.37770414352417, 2.4016101360321045, 0.011450201272964478, -0.6292238831520081];
    end
    [ny, nu] = size(RM);
    nx = nu;
    % do not scale R here because of the observer
    %scale_R = 1e6;
    %RM = RM * scale_R;

    %% Observer
    Lxd_obs = Kfd;
    Lx8_obs = Kfx;    

    S_sp_pinv = pinv([eye(nx)-Ao, -Bo; Co, zeros(ny, nu)]);
    S_sp_pinv = S_sp_pinv(:, nx+1:end);
    S_sp_pinv_x = S_sp_pinv(1:nx,:);
    S_sp_pinv_u = S_sp_pinv(nx+1:end,:);
    
    Apow1 = (Ao(1,1)^1);
    Apow2 = (Ao(1,1)^2);
    Apow3 = (Ao(1,1)^3);
    Apow4 = (Ao(1,1)^4);
    Apow5 = (Ao(1,1)^5);
    Apow6 = (Ao(1,1)^6);
    Apow7 = (Ao(1,1)^7);
    Apow8 = (Ao(1,1)^8);
    
    %% MPC
    horizon = 1;
    u_rate = 0.99*1000; % |(u - awr) * 1e3/1e6| <= 1
    u_max = hardlimits(id_to_cm)*1000; % it's cm_input = |u * 1e3/1e6| <= 5
    J_ = Bo'*P_mpc*Bo+R_mpc;
    L = max(eig(J_)); 
    mu = min(eig(J_));
    J_MPC = eye(size(J_))-J_/L;
    beta_fgm = (sqrt(L) - sqrt(mu)) / (sqrt(L) + sqrt(mu));
    
    q_mat_x0 = Bo'*P_mpc*Ao;
    q_mat_xd = -[Bo'*P_mpc, R_mpc]*[S_sp_pinv_x; S_sp_pinv_u]*(-Cd);
    q_mat_x0 = q_mat_x0/L;
    q_mat_xd = q_mat_xd/L;
    q_mat = [q_mat_x0, q_mat_xd];
    
    %% Rate limiter on VME processors
    a_awr = 2*2*pi;
    g_awr_z = tf([a_awr/(a_awr+2/Ts),a_awr/(a_awr+2/Ts)],...
        [1,(a_awr-2/Ts)/(a_awr+2/Ts)],Ts,'Variable','z^-1');
    sys_awr = eye(nu) * g_awr_z;
    ss_awr = ss(sys_awr);

    %% Get test data
    n_samples = 1000;
    doff = ones(TOT_BPM, n_samples)*10;
    SOFB_setp = 0.9*SOFB_setpoints(id_to_cm)';
    SOFB_setp(SOFB_setp>u_max) = u_max(SOFB_setp>u_max);
    SOFB_setp(SOFB_setp<-u_max) = -u_max(SOFB_setp<-u_max);
    [y_sim,u_sim,obs_y,obs_u,obs_x0,obs_xd,...
     fgm_x0,fgm_xd,fgm_u,fgm_out] = sim_mpc(...
            n_samples, n_delay, doff,...
            Ap, Bp, Cp, ... % Plant
            Ao, Bo, Co, Ad, Cd, Lx8_obs, Lxd_obs,... % Observer
            J_MPC , q_mat, beta_fgm,... % FGM
            u_max , u_rate*ones(nu,1),... % FGM
            id_to_bpm, id_to_cm,...
            ss_awr.A,ss_awr.B,ss_awr.C,...
            SOFB_setp);   
    y_awr = lsim(sys_awr,u_sim(1:n_samples,id_to_cm)/1000,(0:1:n_samples-1)*Ts);
    if true == false
        figure;
        subplot(1,4,1); plot(doff(id_to_bpm,:)'); title('Disturbance [mum]');
        subplot(1,4,2); plot(y_sim(:,id_to_bpm)); title('Output [mum]');
        subplot(1,4,3); plot((u_sim(1:n_samples,id_to_cm)+ones(n_samples,1)*SOFB_setpoints(id_to_cm))/1000); title('Input [A]');
        subplot(1,4,4); plot(u_sim/1000-y_awr); title('AWR [A]');
    end
    if true == true
        figure;
        subplot(1,5,1); plot(doff(id_to_bpm,:)'); title('Disturbance [mum]');
        subplot(1,5,2); plot(y_sim(:,id_to_bpm)); title('Output [mum]');
        subplot(1,5,3); plot((u_sim(1:n_samples,id_to_cm)+ones(n_samples,1)*SOFB_setpoints(id_to_cm))/1000); title('Input with SOFB [A]');
        subplot(1,5,4); plot((u_sim(1:n_samples,id_to_cm))/1000); title('Input [A]');
        subplot(1,5,5); plot(u_sim(1:n_samples,id_to_cm)/1000-y_awr); title('AWR [A]');
    end
    %% FGM Data
    [rows_J_MPC, cols_J_MPC] = size(J_MPC);
    NWORKERS = 6;
    NROWS_W = 32;
    NROWS = NROWS_W * NWORKERS;
    NPAD = NROWS-cols_J_MPC;
    NPAD_x0 = NROWS-nx;
    NPAD_xd = NROWS-ny;
    assert(NPAD>=0);
    assert(NPAD_x0>=0);
    assert(NPAD_xd>=0);

    % zero pad matrix
    J = [J_MPC, zeros(rows_J_MPC, NPAD);...
            zeros(NPAD, NROWS)];
    q_mat_x0_pad = [q_mat_x0, zeros(rows_J_MPC, NPAD_x0);...
                    zeros(NPAD, nx+NPAD_x0)];
    q_mat_xd_pad = [q_mat_xd, zeros(rows_J_MPC, NPAD_xd);...
                    zeros(NPAD, ny+NPAD_xd)];
    q_mat_pad = [q_mat_x0_pad, q_mat_xd_pad];
    
    fid = fopen(sprintf('%s/FGM_data_c6678_NWORKERS6_%s.h',folder_out,lower(dirs{pick_dir})), 'w');
    fprintf(fid, '#ifndef FGM_DATA_C6678_NWORKERS6_%s_H_\n#define FGM_DATA_C6678_NWORKERS6_%s_H_\n\n',...
        upper(dirs{pick_dir}),upper(dirs{pick_dir}));
    fprintf(fid, '#include "mpc/fast_gradient_method.h"\n\n');
    
    fprintf(fid, '#pragma DATA_SECTION(beta_fgm, ".mpc_init");\n');
    fprintf(fid, 'const fgm_float beta_fgm=%.16f;\n', beta_fgm);
    print_dense_C_matrix(fid, [u_max; zeros(NPAD,1)], 'fgm_float', 'ampl_max_vec', true, '.mpc_init', 2);
    print_dense_C_matrix(fid, u_rate*ones(NROWS,1), 'fgm_float', 'rate_max_vec', true, '.mpc_init', 2);
    print_dense_C_matrix(fid, J, 'fgm_float', 'obj_func_matrix', true, '.mpc_init', 2);
    print_dense_C_matrix(fid, q_mat_pad, 'fgm_float', 'obj_func_vector_matrix', true, '.mpc_init', 2);
    %print_dense_C_matrix(fid, [SOFB_setpoints'; zeros(NPAD,1)], 'fgm_float', 'sofb_setp', true, '.mpc_init', 2);
    fprintf(fid, '#endif /* FGM_DATA_C6678_NWORKERS6_%s_H_ */\n',upper(dirs{pick_dir}));
    fclose(fid);
    
    %% FGM test data
    if false
    fgm_x0_pad = [fgm_x0; zeros(NPAD_x0,n_samples)]';
    fgm_xd_pad = [fgm_xd; zeros(NPAD_xd,n_samples)]';
    fgm_u_pad = [fgm_u; zeros(NPAD_x0,n_samples)]';
    fgm_out_pad = [fgm_out; zeros(NPAD_x0,n_samples)]';
        
    fid = fopen(sprintf('%s/FGM_test_data_c6678_NWORKERS6_%s.h',folder_out,lower(dirs{pick_dir})), 'w');
    fprintf(fid, '#ifndef FGM_TEST_DATA_C6678_NWORKERS6_%s_H_\n#define FGM_TEST_DATA_C6678_NWORKERS6_%s_H_\n\n',...
        upper(dirs{pick_dir}),upper(dirs{pick_dir}));
    fprintf(fid, '#include "mpc/fast_gradient_method.h"\n');
    print_dense_C_matrix(fid, fgm_x0_pad, 'fgm_float', 'fgm_x0_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, fgm_xd_pad, 'fgm_float', 'fgm_xd_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, fgm_u_pad, 'fgm_float', 'fgm_u_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, fgm_out_pad, 'fgm_float', 'fgm_out_pad', true, '.mpc_test', 2);
    fprintf(fid, '#endif /* FGM_TEST_DATA_C6678_NWORKERS6_%s_H_ */\n',upper(dirs{pick_dir}));
    fclose(fid);
    end
    
    %% Watchdog data
%     fid = fopen(sprintf('%s/FGM_watchdog_data.h',folder_out), 'w');
%     fprintf(fid, '#ifndef MPC_WATCHDOG_DATA_H_\n#define FGM_WATCHDOG_DATA_H_\n\n');
%     fprintf(fid, '#include "mpc/fast_gradient_method.h"\n\n');
%     print_dense_C_matrix(fid, [u_max; zeros(NPAD,1)], 'fgm_float',...
%         'HARDLIMITS_mA', true, '.fgm_local_data', 2);
%     fprintf(fid, '#endif /* FGM_WATCHDOG_DATA_H_ */\n');
%     fclose(fid);
    
    %% OBS Data
    Lx_pad = [Lx8_obs, zeros(nx, NPAD_xd); zeros(NPAD_x0, NROWS)];
    Ld_pad = [Lxd_obs, zeros(ny, NPAD_xd); zeros(NPAD_xd, NROWS)];
    Cx_pad = [Co, zeros(ny, NPAD_x0); zeros(NPAD_xd, NROWS)];
    
    fid = fopen(sprintf('%s/OBS_data_c6678_NWORKERS6_%s.h',folder_out,lower(dirs{pick_dir})), 'w');
    fprintf(fid, '#ifndef OBS_DATA_C6678_NWORKERS6_%s_H_\n#define OBS_DATA_C6678_NWORKERS6_%s_H_\n\n',...
        upper(dirs{pick_dir}),upper(dirs{pick_dir}));
    fprintf(fid, '#include "mpc/observer.h"\n');
    print_dense_C_matrix(fid, Lx_pad, 'obs_float', 'Lx_pad', true, '.mpc_init', 2);
    print_dense_C_matrix(fid, Ld_pad, 'obs_float', 'Ld_pad', true, '.mpc_init', 2);
    print_dense_C_matrix(fid, Cx_pad, 'obs_float', 'Cx_pad', true, '.mpc_init', 2);
    fprintf(fid, '#endif /* OBS_DATA_C6678_NWORKERS6_%s_H_ */\n',upper(dirs{pick_dir}));
    fclose(fid);    
    
    %% OBS test data
    if false
    obs_u_pad = [obs_u; zeros(NPAD_x0,n_samples)]';
    obs_y_pad = [obs_y; zeros(NPAD_xd,n_samples)]';
    obs_x0_pad = [obs_x0; zeros(NPAD_x0,n_samples)]';
    obs_xd_pad = [obs_xd; zeros(NPAD_xd,n_samples)]';
    
    fid = fopen(sprintf('%s/OBS_test_data_c6678_NWORKERS6_%s.h',folder_out,lower(dirs{pick_dir})), 'w');
    fprintf(fid, '#ifndef OBS_TEST_DATA_C6678_NWORKERS6_%s_H_\n#define OBS_TEST_DATA_C6678_NWORKERS6_%s_H_\n\n',...
        upper(dirs{pick_dir}),upper(dirs{pick_dir}));
    fprintf(fid, '#include "mpc/observer.h"\n');
    print_dense_C_matrix(fid, obs_u_pad, 'obs_float', 'obs_u_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, obs_y_pad, 'obs_float', 'obs_y_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, obs_x0_pad, 'obs_float', 'obs_x0_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, obs_xd_pad, 'obs_float', 'obs_xd_pad', true, '.mpc_test', 2);
    print_dense_C_matrix(fid, SOFB_setpoints, 'obs_float', 'obs_u_pad', true, '.mpc_test', 2);
    fprintf(fid, '#endif /* OBS_TEST_DATA_C6678_NWORKERS6_%s_H_ */\n',upper(dirs{pick_dir}));
    fclose(fid);
    end
    
    %% MPC test data
    network_scaling = -1e3;
    nt = 100;
    ydata = round(y_sim(1:1+nt,:)*1000,0);
    udata = round(u_sim(1+n_delay:1+n_delay+nt,:)*network_scaling,0);
    fid = fopen([folder_out,'MPC_test_data_',lower(dirs{pick_dir}),'.h'], 'w');
    fprintf(fid, sprintf('#ifndef MPC_TEST_DATA_%s_H\n#define MPC_TEST_DATA_%s_H\n',upper(dirs{pick_dir}),upper(dirs{pick_dir})));
    fprintf(fid, '#define MPC_NTEST (%d)\n', nt);
    print_dense_C_matrix(fid, ydata, 'int', 'MPC_TEST_IN', true,'.mpc_unit_test', 2);
    print_dense_C_matrix(fid, udata, 'int', 'MPC_TEST_OUT', true,'.mpc_unit_test', 2);
    fprintf(fid, '#endif\n');
    fclose(fid);

    %% Storage Ring Config
    fid = fopen(sprintf('%s/MPC_storage_ring_config_%s.h',folder_out,lower(dirs{pick_dir})), 'w');
    fprintf(fid, '#ifndef MPC_STORAGE_RING_CONFIG_%s_H_\n#define MPC_STORAGE_RING_CONFIG_%s_H_\n\n',...
        upper(dirs{pick_dir}),upper(dirs{pick_dir}));
    fprintf(fid, '#define MPC_NU_ (%d)\n', nu);
    fprintf(fid, '#define MPC_NY_ (%d)\n\n', ny);
    print_dense_C_matrix(fid, id_to_bpm-1, 'int', 'MPC_ID_TO_BPM', true, '.mpc_shared', 2);
    print_dense_C_matrix(fid, id_to_cm-1, 'int', 'MPC_CM_TO_ID', true, '.mpc_shared', 2);
    fprintf(fid, '#endif /* MPC_STORAGE_RING_CONFIG_%s_H_ */\n',upper(dirs{pick_dir}));
    fclose(fid);
    
    %% AWR
    [num_awr, den_awr] = tfdata(g_awr_z);
    generate_filter_from_TF(folder_out, NROWS_W, num_awr{:}, den_awr{:}, 'awr',...
        'float', 0, '.mpc_awr', 64, 1e20);

end