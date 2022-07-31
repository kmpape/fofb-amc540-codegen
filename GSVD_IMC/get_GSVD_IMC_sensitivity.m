function [Sminmax_x, Sminmax_y] = get_GSVD_IMC_sensitivity(RMorigx, RMorigy, bws, bwf, w_Hz)
addpath('..')
if ~exist('bws','var')
    bws = 100*2*pi;
end
if ~exist('bwf','var')
    bwf = 200*2*pi;
end

%% Configure Diamond-I Storage Ring
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
n_delay = 8; % number of delay time steps [-]
z = tf('z');
s = tf('s');
aIx = 2*pi*500;
aIy = 2*pi*700;
gx = aIx/(s+aIx);
gy = aIy/(s+aIy);
T_delay = n_delay*Ts; % delay [s]
tf_delay = exp(-s*T_delay);

%% Mid-Ranging IMC
bw_allx = bwf; % overall desired bandwidth [rad/s]
bw_sx = bws; % slow actuators desired bandwidth [rad/s]
T_tiso_mpx = bw_allx/(s+bw_allx);
T_siso_mpx = bw_sx/(s+bw_sx);

bw_ally = bwf; % overall desired bandwidth [rad/s]
bw_sy = bws; % slow actuators desired bandwidth [rad/s]
T_tiso_mpy = bw_ally/(s+bw_ally);
T_siso_mpy = bw_sy/(s+bw_sy);

T_tisox = T_tiso_mpx * tf_delay;
T_sisox = T_siso_mpx * tf_delay;
T_tisoy = T_tiso_mpy * tf_delay;
T_sisoy = T_siso_mpy * tf_delay;

qsx = T_siso_mpx / gx;
qfx = (T_tiso_mpx - T_siso_mpx) / gx;
qsy = T_siso_mpy / gy;
qfy = (T_tiso_mpy - T_siso_mpy) / gy;

%% GSVD
[Usx,Ufx,Xx,C,S] = gsvd(Rsx', Rfx');
Ssx = C';
Sfx = S';
Ffx = Xx*pinv(Xx*blkdiag(eye(rank(Rfx)), zeros(ns-rank(Rfx))));  
mu = 1;
Gx = Xx*((Xx'*Xx+mu*eye(ny))\Xx');
Ksx = -(Usx*inv(Ssx)/Xx)*Gx;
Kfx = -((Ufx*pinv(Sfx)/Xx)*Ffx)*Gx;
Psx = -Gx\Rsx;
Pfx = -Gx\Rfx;

[Usy,Ufy,Xy,C,S] = gsvd(Rsy', Rfy');
Ssy = C';
Sfy = S';
Ffy = Xy*pinv(Xy*blkdiag(eye(rank(Rfy)), zeros(ns-rank(Rfy)))); 
mu = 1;
Gy = Xy*((Xy'*Xy+mu*eye(ny))\Xy');
Ksy = -(Usy*inv(Ssy)/Xy)*Gy;
Kfy = -((Ufy*pinv(Sfy)/Xy)*Ffy)*Gy;
Psy = -Gy\Rsy;
Pfy = -Gy\Rfy;

%%
nw = length(w_Hz);

T_tisox_w = freqresp(T_tisox, w_Hz, 'Hz'); T_tisox_w = T_tisox_w(:);
T_sisox_w = freqresp(T_sisox, w_Hz, 'Hz'); T_sisox_w = T_sisox_w(:);
T_tisoy_w = freqresp(T_tisoy, w_Hz, 'Hz'); T_tisoy_w = T_tisoy_w(:);
T_sisoy_w = freqresp(T_sisoy, w_Hz, 'Hz'); T_sisoy_w = T_sisoy_w(:);

qsx_w = freqresp(qsx, w_Hz, 'Hz'); qsx_w = qsx_w(:);
qfx_w = freqresp(qfx, w_Hz, 'Hz'); qfx_w = qfx_w(:);
gsx_w = freqresp(gx, w_Hz, 'Hz'); gsx_w = gsx_w(:);
gfx_w = freqresp(gx, w_Hz, 'Hz'); gfx_w = gfx_w(:);

qsy_w = freqresp(qsy, w_Hz, 'Hz'); qsy_w = qsy_w(:);
qfy_w = freqresp(qfy, w_Hz, 'Hz'); qfy_w = qfy_w(:);
gsy_w = freqresp(gy, w_Hz, 'Hz'); gsy_w = gsy_w(:);
gfy_w = freqresp(gy, w_Hz, 'Hz'); gfy_w = gfy_w(:);


tmp_Qsx = Usx*inv(Ssx)/Xx;
tmp_Qfx = Ufx*pinv(Sfx)/Xx;
tmp_Qsy = Usy*inv(Ssy)/Xy;
tmp_Qfy = Ufy*pinv(Sfy)/Xy;
Sminmax_x = zeros(2,nw);
Sminmax_y = zeros(2,nw);
I = eye(ny);
for i = 1 : nw
    
    Qsx = tmp_Qsx*qsx_w(i);
    Qfx = tmp_Qfx*qfx_w(i);
    Qsy = tmp_Qsy*qsy_w(i);
    Qfy = tmp_Qfy*qfy_w(i);
    
    Psx = Rsx*gsx_w(i);
    Pfx = Rfx*gfx_w(i);
    Psy = Rsy*gsy_w(i);
    Pfy = Rfy*gfy_w(i);
    
    S_matx = I-(Psx*Qsx+Pfx*Qfx*Ffx)*inv(I+(Gx-I)*(Psx*Qsx+Pfx*Qfx*Ffx))*Gx;
    S_maty = I-(Psy*Qsy+Pfy*Qfy*Ffy)*inv(I+(Gy-I)*(Psy*Qsy+Pfy*Qfy*Ffy))*Gy;
    
    Sminmax_x(1,i) = min(svd(S_matx));
    Sminmax_x(2,i) = max(svd(S_matx));
    Sminmax_y(1,i) = min(svd(S_maty));
    Sminmax_y(2,i) = max(svd(S_maty));
end

end



