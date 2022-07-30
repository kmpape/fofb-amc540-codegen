function [n_delay, Ax, Bx, Cx, Dx, Ay, By, Cy, Dy,...
          Ksx, Acsx, Bcsx, Ccsx, Dcsx, Ksy, Acsy, Bcsy, Ccsy, Dcsy,...
          Kfx, Acfx, Bcfx, Ccfx, Dcfx, Kfy, Acfy, Bcfy, Ccfy, Dcfy,...
          Psx, Amsx, Bmsx, Cmsx, Dmsx, Psy, Amsy, Bmsy, Cmsy, Dmsy,...
          Pfx, Amfx, Bmfx, Cmfx, Dmfx, Pfy, Amfy, Bmfy, Cmfy, Dmfy,...
          id_to_bpm, slow_to_id, fast_to_id, network_scaling] = get_GSVD_IMC_controller(RMorigx, RMorigy, bws, bwf)
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
tf_DIx = aIx/(s+aIx);
tf_DIy = aIy/(s+aIy);
gI_mp_zx = c2d(tf_DIx, Ts, 'zoh');
gI_mp_zy = c2d(tf_DIy, Ts, 'zoh');

% input is negated on PMCs, so multiply with -1 here
minus_one = -1;

%% Mid-Ranging IMC
bw_allx = bwf; % overall desired bandwidth [rad/s]
bw_sx = bws; % slow actuators desired bandwidth [rad/s]
T_tiso_mpx = bw_allx/(s+bw_allx);
T_siso_mpx = bw_sx/(s+bw_sx);

bw_ally = bwf; % overall desired bandwidth [rad/s]
bw_sy = bws; % slow actuators desired bandwidth [rad/s]
T_tiso_mpy = bw_ally/(s+bw_ally);
T_siso_mpy = bw_sy/(s+bw_sy);

qs_zx = c2d(T_siso_mpx / tf_DIx, Ts, 'zoh');
qf_zx = c2d((T_tiso_mpx - T_siso_mpx) / tf_DIx, Ts, 'zoh');
qs_zy = c2d(T_siso_mpy / tf_DIy, Ts, 'zoh');
qf_zy = c2d((T_tiso_mpy - T_siso_mpy) / tf_DIy, Ts, 'zoh');

%% GSVD
[Us,Uf,X,C,S] = gsvd(Rsx', Rfx');
Ss = C';
Sf = S';
Ff = X*pinv(X*blkdiag(eye(rank(Rfx)), zeros(ns-rank(Rfx))));  
mu = 1;
G = X*((X'*X+mu*eye(ny))\X');
Ksx = -(Us*inv(Ss)/X)*G;
Kfx = -((Uf*pinv(Sf)/X)*Ff)*G;
Psx = -G\Rsx;
Pfx = -G\Rfx;

[Us,Uf,X,C,S] = gsvd(Rsy', Rfy');
Ss = C';
Sf = S';
Ff = X*pinv(X*blkdiag(eye(rank(Rfy)), zeros(ns-rank(Rfy)))); 
mu = 1;
G = X*((X'*X+mu*eye(ny))\X');
Ksy = -(Us*inv(Ss)/X)*G;
Kfy = -((Uf*pinv(Sf)/X)*Ff)*G;
Psy = -G\Rsy;
Pfy = -G\Rfy;

%%
[Ax, Bx, Cx ,Dx] = ssdata(RMorigx .* gI_mp_zx);
[Acsx,Bcsx,Ccsx,Dcsx] = ssdata(eye(ns) .* c2d(T_siso_mpx / tf_DIx, Ts, 'zoh'));
[Acfx,Bcfx,Ccfx,Dcfx] = ssdata(eye(nf) .* c2d((T_tiso_mpx - T_siso_mpx) / tf_DIx, Ts, 'zoh'));
[Amsx,Bmsx,Cmsx,Dmsx] = ssdata(eye(ny) .* gI_mp_zx .* z^(-n_delay));
[Amfx,Bmfx,Cmfx,Dmfx] = ssdata(eye(ny) .* gI_mp_zx .* z^(-n_delay));

[Ay, By, Cy ,Dy] = ssdata(RMorigy .* gI_mp_zy);
[Acsy,Bcsy,Ccsy,Dcsy] = ssdata(eye(ns) .* c2d(T_siso_mpy / tf_DIy, Ts, 'zoh'));
[Acfy,Bcfy,Ccfy,Dcfy] = ssdata(eye(nf) .* c2d((T_tiso_mpy - T_siso_mpy) / tf_DIy, Ts, 'zoh'));
[Amsy,Bmsy,Cmsy,Dmsy] = ssdata(eye(ny) .* gI_mp_zy .* z^(-n_delay));
[Amfy,Bmfy,Cmfy,Dmfy] = ssdata(eye(ny) .* gI_mp_zy .* z^(-n_delay));

network_scaling = 1e6;
Ksx = Ksx / network_scaling;
Kfx = Kfx / network_scaling;
Psx = Psx * network_scaling;
Pfx = Pfx * network_scaling;

Ksy = Ksy / network_scaling;
Kfy = Kfy / network_scaling;
Psy = Psy * network_scaling;
Pfy = Pfy * network_scaling;

end



