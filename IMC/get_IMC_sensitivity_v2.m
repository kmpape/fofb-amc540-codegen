function [Sminmax_x, Sminmax_y] = get_IMC_sensitivity_v2(RMorigx, RMorigy, bw, w_Hz, config_num, n_delay)

addpath('..')
plot_bode = false;

%% Configure Diamond-I Storage Ring

[ny_x,nu_x] = size(RMorigx);
[ny_y,nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
if config_num == 5 % use same storage ring config as GSVD-IMC
    [id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v5(RMorigx,RMorigy,false);
else % full storage ring
    assert(0);
end

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
Ix = eye(ny_x);
Iy = eye(ny_y);
for i = 1 : nw
    
    Cx = Kx*cx_w(i);
    Cy = Ky*cy_w(i);
    
    Px = RMx*gx_w(i);
    Py = RMy*gy_w(i);
    
    S_matx = inv(Ix+Px*Cx);
    S_maty = inv(Iy+Py*Cy);
    
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
