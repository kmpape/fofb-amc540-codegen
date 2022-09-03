function [network_scaling,...
          Acx, Bcx, Ccx, Dcx, Ax, Bx, Cx, Dx,...
          Acy, Bcy, Ccy, Dcy, Ay, By, Cy, Dy] = get_IMC_controller_cfgv4_square(RMorigx, RMorigy,...
                                                                                bw, n_delay)
% Returns the controller in standard FB form with the corresponding SS
% system
addpath('..')

%% Configure Diamond-I Storage Ring

[ny_x,nu_x] = size(RMorigx);
[ny_y,nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
[id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v4(RMorigx,RMorigy,true);

RMx = RMorigx(id_to_bpm_x, id_to_cm_x);
RMy = RMorigy(id_to_bpm_y, id_to_cm_y);

ny_x = length(id_to_bpm_x);
ny_y = length(id_to_bpm_y);

%% Actuators
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time[s]
z = tf('z');
s = tf('s');
aIx = 2*pi*500;
aIy = 2*pi*700;
tf_DIx = aIx/(s+aIx);
tf_DIy = aIy/(s+aIy);
gI_mp_zx = c2d(tf_DIx, Ts, 'zoh');
gI_mp_zy = c2d(tf_DIy, Ts, 'zoh');

%% IMC
[U, S, V] = svd(RMx, 'econ');
MU = 1.0*eye(ny_x);
E = S / (S.^2+MU);
Kx = V*E*U';

[U, S, V] = svd(RMy, 'econ');
MU = 1.0*eye(ny_y);
E = S / (S.^2+MU);
Ky = V*E*U';

z = tf('z', Ts);
s = tf('s');
%bw = 2*pi*200;
abw = exp(-bw*Ts);
T_mp_z = (1-abw) / (1-z^(-1)*abw) * z^(-1);
q_zx = T_mp_z / gI_mp_zx;
q_zy = T_mp_z / gI_mp_zy;

c_zx = q_zx / (1 - T_mp_z*z^(-n_delay));
c_zy = q_zy / (1 - T_mp_z*z^(-n_delay));


%%
network_scaling = 1e6;
[Ax, Bx, Cx ,Dx] = ssdata(RMorigx .* gI_mp_zx);
[Ay, By, Cy ,Dy] = ssdata(RMorigy .* gI_mp_zy);
[Acx, Bcx, Ccx, Dcx] = ssdata(-(Kx/network_scaling) .* c_zx);
[Acy, Bcy, Ccy, Dcy] = ssdata(-(Ky/network_scaling) .* c_zy);


end
