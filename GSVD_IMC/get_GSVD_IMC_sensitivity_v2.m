function [Sminmax_x, Sminmax_y] = get_GSVD_IMC_sensitivity_v2(RMorigx, RMorigy, bwsx, bwfx, bwsy, bwfy, n_delay, w_Hz, plot_bode, n_delay_plant, mu)
addpath('..')
if ~exist('n_delay','var') || isempty(n_delay)
    n_delay = 8;
end
if ~exist('n_delay_plant','var') || isempty(n_delay_plant)
    n_delay_plant = n_delay;
end
if ~exist('plot_bode','var') || isempty(plot_bode)
    plot_bode = false;
end
if ~exist('mu','var') || isempty(mu)
    mu = 1;
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

id_to_cm = sort(union(slow_to_id,fast_to_id), 'asc');
Rx = RMorigx(id_to_bpm, id_to_cm);
Ry = RMorigy(id_to_bpm, id_to_cm);

%% Actuators
Fs = 10*10^3; % sample frequency [Hz]
Ts = 1/Fs; % sample time[s]
z = tf('z');
s = tf('s');
aIx = 2*pi*500;
aIy = 2*pi*700;
gx_mp = aIx/(s+aIx);
gy_mp = aIy/(s+aIy);
T_delay = n_delay*Ts; % delay [s]
T_delay_plant = n_delay_plant*Ts; % delay [s]
tf_delay = exp(-s*T_delay);
tf_delay_plant = exp(-s*T_delay_plant);
gx = gx_mp * tf_delay;
gy = gy_mp * tf_delay;
gx_plant = gx_mp * tf_delay_plant;
gy_plant = gy_mp * tf_delay_plant;

bode_opt = bodeoptions();
bode_opt.FreqUnits = 'Hz';
bode_opt.XLim = [0.01, 10^4];
if false
    figure;
    bodeplot(gx, gy, bode_opt);
    axes_handles = findall(gcf, 'type', 'axes');
    legend(axes_handles(3),'X','Y','Location', 'NorthEast')
    title('Actuators');
end

%% Mid-Ranging IMC
bw_allx = bwfx; % overall desired bandwidth [rad/s]
bw_sx = bwsx; % slow actuators desired bandwidth [rad/s]
T_tiso_mpx = bw_allx/(s+bw_allx);
T_siso_mpx = bw_sx/(s+bw_sx);

bw_ally = bwfy; % overall desired bandwidth [rad/s]
bw_sy = bwsy; % slow actuators desired bandwidth [rad/s]
T_tiso_mpy = bw_ally/(s+bw_ally);
T_siso_mpy = bw_sy/(s+bw_sy);

T_tisox = T_tiso_mpx * tf_delay;
T_sisox = T_siso_mpx * tf_delay;
T_tisoy = T_tiso_mpy * tf_delay;
T_sisoy = T_siso_mpy * tf_delay;

S_tisox = 1 - T_tisox;
S_sisox = 1 - T_sisox;
S_tisoy = 1 - T_tisoy;
S_sisoy = 1 - T_sisoy;

qsx = T_siso_mpx / gx_mp;
qfx = (T_tiso_mpx - T_siso_mpx) / gx_mp;
qsy = T_siso_mpy / gy_mp;
qfy = (T_tiso_mpy - T_siso_mpy) / gy_mp;

if false
    figure;
    bodeplot(T_tisox, T_sisox, bode_opt);
    axes_handles = findall(gcf, 'type', 'axes');
    legend(axes_handles(3),'TISO X \& Y','SISO X \& Y','Location', 'SouthWest');
    title('Complementary Sensitivity');
            
    figure;
    bodeplot(S_tisox, S_sisox, bode_opt);
    axes_handles = findall(gcf, 'type', 'axes');
    legend(axes_handles(3),'TISO X \& Y','SISO X \& Y','Location', 'SouthWest');
    title('Output Sensitivity');
end
    
if false
    figure;
    bodeplot(qsx, qfx, qsy, qfy, bode_opt);
    axes_handles = findall(gcf, 'type', 'axes');
    legend(axes_handles(3),'$q_{s,x}$', '$q_{f,x}$', '$q_{s,y}$', '$q_{f,y}$','Location', 'SouthWest');
    title('IMC Controller');
end

%% GSVD
%mu = 1;

[Usx,Ufx,Xx,C,S] = gsvd(Rsx', Rfx');
% [UX,~,~] = svd(Xx);
% m_crit = size(Xx,1);
% w_vec = [ones(m_crit,1); ones(size(Xx,1)-m_crit,1)/100];
% W = UX * diag(w_vec) * UX' / sum(w_vec) * size(Xx,1);
Ssx = C';
Sfx = S';
Ffx = Xx*pinv(Xx*blkdiag(eye(rank(Rfx)), zeros(ns-rank(Rfx))));  
Gx = Xx*((Xx'*Xx+mu*eye(ny))\Xx');
%Gx = Xx*((Xx'*W*Xx+mu*eye(ny))\Xx'*W');
% Ksx = -(Usx*inv(Ssx)/Xx)*Gx;
% Kfx = -((Ufx*pinv(Sfx)/Xx)*Ffx)*Gx;
% Psx = -Gx\Rsx;
% Pfx = -Gx\Rfx;

[Usy,Ufy,Xy,C,S] = gsvd(Rsy', Rfy');
Ssy = C';
Sfy = S';
Ffy = Xy*pinv(Xy*blkdiag(eye(rank(Rfy)), zeros(ns-rank(Rfy))));

Gy = Xy*((Xy'*Xy+mu*eye(ny))\Xy');
% Ksy = -(Usy*inv(Ssy)/Xy)*Gy;
% Kfy = -((Ufy*pinv(Sfy)/Xy)*Ffy)*Gy;
% Psy = -Gy\Rsy;
% Pfy = -Gy\Rfy;

%%
nw = length(w_Hz);

T_tisox_w = freqresp(T_tisox, w_Hz, 'Hz'); T_tisox_w = T_tisox_w(:);
T_sisox_w = freqresp(T_sisox, w_Hz, 'Hz'); T_sisox_w = T_sisox_w(:);
T_tisoy_w = freqresp(T_tisoy, w_Hz, 'Hz'); T_tisoy_w = T_tisoy_w(:);
T_sisoy_w = freqresp(T_sisoy, w_Hz, 'Hz'); T_sisoy_w = T_sisoy_w(:);

S_tisox_w = freqresp(S_tisox, w_Hz, 'Hz'); S_tisox_w = S_tisox_w(:);
S_sisox_w = freqresp(S_sisox, w_Hz, 'Hz'); S_sisox_w = S_sisox_w(:);
S_tisoy_w = freqresp(S_tisoy, w_Hz, 'Hz'); S_tisoy_w = S_tisoy_w(:);
S_sisoy_w = freqresp(S_sisoy, w_Hz, 'Hz'); S_sisoy_w = S_sisoy_w(:);

qsx_w = freqresp(qsx, w_Hz, 'Hz'); qsx_w = qsx_w(:);
qfx_w = freqresp(qfx, w_Hz, 'Hz'); qfx_w = qfx_w(:);
gsx_w = freqresp(gx_plant, w_Hz, 'Hz'); gsx_w = gsx_w(:);
gfx_w = freqresp(gx_plant, w_Hz, 'Hz'); gfx_w = gfx_w(:);

qsy_w = freqresp(qsy, w_Hz, 'Hz'); qsy_w = qsy_w(:);
qfy_w = freqresp(qfy, w_Hz, 'Hz'); qfy_w = qfy_w(:);
gsy_w = freqresp(gy_plant, w_Hz, 'Hz'); gsy_w = gsy_w(:);
gfy_w = freqresp(gy_plant, w_Hz, 'Hz'); gfy_w = gfy_w(:);


tmp_Qsx = Usx*inv(Ssx)/Xx;
tmp_Qfx = Ufx*pinv(Sfx)/Xx;
tmp_Qsy = Usy*inv(Ssy)/Xy;
tmp_Qfy = Ufy*pinv(Sfy)/Xy;
Sminmax_x = zeros(2,nw);
Sminmax_y = zeros(2,nw);

if true
[URx,~,~] = svd(Rx);
[URy,~,~] = svd(Ry);
Smode_x = zeros(ny,nw);
Smode_y = zeros(ny,nw);
end

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
    
    if true
        Smode_x(:,i) = vecnorm(URx'*S_matx*URx, 2, 2);
        Smode_y(:,i) = vecnorm(URy'*S_maty*URy, 2, 2);
    end
end

if plot_bode
    figure;
    subplot(1,2,1);
    semilogx(w_Hz, 20*log10(Sminmax_x(1,:)), 'r', w_Hz, 20*log10(Sminmax_x(2,:)), 'r', ...
        w_Hz, 20*log10(abs(S_sisox_w)), 'b--', w_Hz, 20*log10(abs(S_tisox_w)), 'b--');
    title('X')
    subplot(1,2,2);
    semilogx(w_Hz, 20*log10(Sminmax_y(1,:)), 'r', w_Hz, 20*log10(Sminmax_y(2,:)), 'r', ...
        w_Hz, 20*log10(abs(S_sisoy_w)), 'b--', w_Hz, 20*log10(abs(S_tisoy_w)), 'b--');
    title('Y')    
end

if plot_bode
    figure;
    subplot(1,2,1);
    h=surf(w_Hz, 1:1:ny, 20*log10(Smode_x));
    view(0,90); set(h,'LineStyle','none'); set(gca,'xscale','log')
    caxis([-80,5]);
    xlim([w_Hz(1) w_Hz(end)]);
    ylim([1 ny]);
    colorbar; title('Mode Space - X');
    subplot(1,2,2);
    h=surf(w_Hz, 1:1:ny, 20*log10(Smode_y));
    view(0,90); set(h,'LineStyle','none'); set(gca,'XScale','log');
    caxis([-80,5]);
    xlim([w_Hz(1) w_Hz(end)]);
    ylim([1 ny]);
    colorbar; title('Mode Space - Y');
end

end



