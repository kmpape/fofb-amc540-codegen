clear all
close all
clc

%%
folder = 'myfolder';
writepath = sprintf('/path_to_folder/%s/', folder); 
if exist(writepath,'dir') ~= 7
    mkdir(writepath);
end

BPM_disable = [];
id_cm = 1:172; % find the magnets on these IDs when BPMs are renumbered
id_bpm_ = 257:429; % find the BPMs on these IDs 
id_bpm = setdiff(id_bpm_,256+BPM_disable); 
ind_cm = 1:172;
ind_bpm = 173:173+length(id_bpm)-1;


%% TIMES
% Format:
% t_vec = {tcell1, tcell2,...}
% tcelli = {start_time,end_time, dt, comment}
% dt: Data from start_time to end_time is split into bits of length dt.
%     Every bit of length dt is saved into a separate file 


% Use this array for the case that BPMs are renumbered
dt = 10;
t_vec = {...
{datetime(2022,9,28,00,25,00), datetime(2022,9,28,00,26,30),dt, '300 mA FOFB off'},... %1
{datetime(2022,9,28,00,30,00), datetime(2022,9,28,00,50,00),dt, '300 mA MPC v1.0.3-xi steady-state'},... %2
{datetime(2022,9,28,01,05,00), datetime(2022,9,28,01,25,00),dt, '300 mA IMC v1.0.4-delta steady-state'},... %3
{datetime(2022,9,28,01,40,00), datetime(2022,9,28,02,00,00),dt, '300 mA GSVD v1.0.4-gamma steady-state'},... %4
};

% Use this array for the case that BPMs are not renumbered
t_vec_standard = {...
{datetime(2022,9,28,00,24,00), datetime(2022,9,28,00,24,40),dt, '300 mA FOFB on'},... %
};

%% USE THIS FOR THE CASE THAT BPMs ARE RENUMBERED
write_mat = true;
dirs = {'X','Y'};

for j = 1 : length(t_vec)
    fprintf("At j=%d\n",j);
    
    tend = t_vec{j}{2};    
    t1 = t_vec{j}{1};
    t2 = t1+seconds(t_vec{j}{3});
    while (t2 <= tend)
        t1_tmp = t1;
        t2_tmp = t2;
        t1_tmp.TimeZone = 'Europe/London';
        t2_tmp.TimeZone = 'Europe/London';
        data_ = fa_load([datenum(t1) datenum(t2)], [id_cm id_bpm], 'F');
        data_pmc = fa_load([datenum(t1) datenum(t2)], [201], 'F');
        [sofb_x, sofb_y] = get_SOFB(t1_tmp);
        [gap_data,gap_IDs] = get_GAP(t1_tmp, t2_tmp);
        t1_str = datestr(t1,'ddmmyyyy_HHMMSS');
    
        n_samples = size(data_.data,3);
        data = data_.data(:,:,1:n_samples);
        decim = 1;
        n = floor(n_samples/decim);
        F_S = 10072;
        T_S = 1/F_S;
        t = (0:1:n_samples-1)*T_S;
        tn = t(1:decim:n_samples);

        for p = 1:2
            if p == 1
                sofb = sofb_x;
            else
                sofb = sofb_y;
            end
            y_pmc = squeeze(data_pmc.data(2,1,:));
            y = squeeze(data(p,ind_bpm,:))*1e-3; % micrometers
            u = -squeeze(data(p,ind_cm,:))/1e6; % Amperes

            if write_mat == true
                fname = sprintf('%s_data_%s.mat',t1_str,dirs{p});
                tmstmp = data_.timestamp;
                tt_vec = t_vec{j};
                save([writepath,fname], 'y', 'u', 't', 'y_pmc','tmstmp', 'tt_vec', 'sofb', 'gap_data', 'gap_IDs');
            end
        end
        t1 = t2;
        t2 = t2 + seconds(t_vec{j}{3});
    end
end

%% USE THIS FOR THE CASE THAT BPMs ARE NOT RENUMBERED
write_mat = true;
dirs = {'X','Y'};
BPM_disable = [];
id_bpm = 1:173; % overwrite this
ind_bpm = id_bpm;
for j = 1 : length(t_vec_standard)
    fprintf("At SFB j=%d\n",j);
    
    tend = t_vec_standard{j}{2};    
    t1 = t_vec_standard{j}{1};
    t2 = t1+seconds(t_vec_standard{j}{3});
    while (t2 <= tend)
    
        data_ = fa_load([datenum(t1) datenum(t2)], [id_bpm], 'F');
        data_pmc = fa_load([datenum(t1) datenum(t2)], [201], 'F');
        [sofb_x, sofb_y] = get_SOFB(t1);
        t1_str = datestr(t1,'ddmmyyyy_HHMMSS');
    
        n_samples = size(data_.data,3);
        data = data_.data(:,:,1:n_samples);
        decim = 1;
        n = floor(n_samples/decim);
        F_S = 10072;
        T_S = 1/F_S;
        t = (0:1:n_samples-1)*T_S;
        tn = t(1:decim:n_samples);

        for p = 1:2
            if p == 1
                sofb = sofb_x;
            else
                sofb = sofb_y;
            end
            y_pmc = squeeze(data_pmc.data(2,1,:));
            y = squeeze(data(p,ind_bpm,:))*1e-3; % micrometers

            if write_mat == true
                fname = sprintf('%s_data_sfb_%s.mat',t1_str,dirs{p});
                tmstmp = data_.timestamp;
                tt_vec = t_vec_standard{j};
                save([writepath,fname], 'y', 't', 'y_pmc','tmstmp', 'tt_vec', 'sofb');
            end
        end
        t1 = t2;
        t2 = t2 + seconds(t_vec_standard{j}{3});
    end
end