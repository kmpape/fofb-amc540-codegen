function [Sminmax_x, Sminmax_y] = get_IMC_sensitivity(RMorigx, RMorigy, bw, w_Hz, full_config, n_delay)

addpath('..')
if ~exist('full_config','var')
    full_config = false;
end
if ~exist('n_delay','var')
    n_delay = 8;
end
plot_bode = false;

%% Configure Diamond-I Storage Ring

[ny_x,nu_x] = size(RMorigx);
[ny_y,nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
if (~full_config) % use same storage ring config as GSVD-IMC
    [id_to_bpm, slow_to_id, ~] = diamond_I_configuration_v3(RMorigx,RMorigy);
else % full storage ring
    id_to_bpm = 1:1:173;
    bad_bpm = [76,79];
    id_to_bpm(bad_bpm) = [];
    slow_to_id = 1:1:172;
end
id_to_cm = slow_to_id;

RMx = RMorigx(id_to_bpm, slow_to_id);
RMy = RMorigy(id_to_bpm, slow_to_id);

ny = length(id_to_bpm);

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
MU = 1.0*eye(ny);
E = S / (S.^2+MU);
Kx = V*E*U';

[U, S, V] = svd(RMy, 'econ');
MU = 1.0*eye(ny);
E = S / (S.^2+MU);
Ky = V*E*U';

z = tf('z', Ts);
s = tf('s');
%bw = 2*pi*200;
abw = exp(-bw*Ts);
T_mp_z = (1-abw) / (z-abw);
q_zx = T_mp_z / gI_mp_zx;
q_zy = T_mp_z / gI_mp_zy;

c_zx = q_zx / (1 - T_mp_z*z^(-n_delay));
c_zy = q_zy / (1 - T_mp_z*z^(-n_delay));

if plot_bode
    S = 1 - T_mp_z*z^(-n_delay);
    S2 = 1/(1+ c_zx*gI_mp_zx*z^(-n_delay));
    bode_opt = bodeoptions();
    bode_opt.FreqUnits = 'Hz';
    bode_opt.XLim = [0.1, 10^4];
    figure;
    bodeplot(S, S2, bode_opt);
    axes_handles = findall(gcf, 'type', 'axes');
    legend(axes_handles(3),'S','S2','Location', 'NorthEast')
    title('Sensitivity');
end


%%
nw = length(w_Hz);
cx_w = freqresp(c_zx, w_Hz, 'Hz'); cx_w = cx_w(:);
gx_w = freqresp(gI_mp_zx*z^(-n_delay), w_Hz, 'Hz'); gx_w = gx_w(:);

cy_w = freqresp(c_zy, w_Hz, 'Hz'); cy_w = cy_w(:);
gy_w = freqresp(gI_mp_zy*z^(-n_delay), w_Hz, 'Hz'); gy_w = gy_w(:);

Sminmax_x = zeros(2,nw);
Sminmax_y = zeros(2,nw);
I = eye(ny);
for i = 1 : nw
    
    Cx = Kx*cx_w(i);
    Cy = Ky*cy_w(i);
    
    Px = RMx*gx_w(i);
    Py = RMy*gy_w(i);
    
    S_matx = inv(I+Px*Cx);
    S_maty = inv(I+Py*Cy);
    
    Sminmax_x(1,i) = min(svd(S_matx));
    Sminmax_x(2,i) = max(svd(S_matx));
    Sminmax_y(1,i) = min(svd(S_maty));
    Sminmax_y(2,i) = max(svd(S_maty));
end

if plot_bode
    Sw = freqresp(S, w_Hz, 'Hz'); Sw = Sw(:);
    figure;
    semilogx(w_Hz, 20*log10(Sminmax_x(1,:)), 'b', w_Hz, 20*log10(Sminmax_x(2,:)), 'b',...
        w_Hz, 20*log10(Sminmax_y(1,:)), 'r', w_Hz, 20*log10(Sminmax_y(1,:)), 'r',...
        w_Hz, 20*log10(abs(Sw)), 'g--');
end

end
