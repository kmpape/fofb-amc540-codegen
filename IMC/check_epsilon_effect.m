

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Model = First Order plus delay
% (1) Vertical 
Ts= 1e-4;
tdy = 700e-6; 
dy = ceil(tdy/Ts);
taudy = tdy - (dy-1)*Ts;
zetay = 1/(1*tdy); 
ly = exp(-zetay*Ts);

% Different dynamics 
% ay(1) = Regular correctors
% ay(2) = Cell2 fast correctors
% ay(3) = Cell2 slow correctors
ay = [2*pi*700; 2*pi*2000; 2*pi*10];
for ii = 1:length(ay)
      ay1(ii) = exp(-ay(ii)*Ts);
      by0(ii) = 1-exp(-ay(ii)*(Ts-taudy)); % why not: 1-exp(-ay(ii)*(Ts-taudy))? -> ~=0...
      by1(ii) = exp(-ay(ii)*(Ts-taudy))-ay1(ii);
      gzy(ii) = tf([zeros(1,dy) by0(ii) by1(ii)], [1 -ay1(ii)], Ts, 'variable', 'z^-1');
end


% (2) Horizontal 
Ts= 1e-4;
tdx = 700e-6; 
dx = ceil(tdx/Ts);
taudx = tdx - (dx-1)*Ts;
zetax = 1/(1*tdx); 
lx = exp(-zetax*Ts);


% Different dynamics 
% ax(1) = Regular correctors
% ax(2) = Cell2 fast correctors
% ax(3) = Cell2 slow correctors
ax = [2*pi*500; 2*pi*2000; 2*pi*20];
for ii = 1:length(ax)
      ax1(ii) = exp(-ax(ii)*Ts);
      bx0(ii) = 0;
      bx1(ii) = exp(-ax(ii)*(Ts-taudx))-ax1(ii);
      gzx(ii) =tf([zeros(1,dx) bx0(ii) bx1(ii)],[1 -ax1(ii)],Ts,'variable','z^-1');
end

%% DYNAMIC CONTROLLER
ep = 1e-3;
%Horizontal SISO Controller
qzx = minreal(tf([1-lx],[1 -lx],Ts,'variable','z^-1')*tf([1 -ax1(1)],[bx0(1)+bx1(1)],Ts,'variable','z^-1'));
epi = 1/(1-ep);
czx = ((1-lx)/(bx0(1)+bx1(1)))*tf([1 -ax1(1)],[epi -lx*epi 0 0 0 0 0 0 -(1-lx)],Ts,'variable','z^-1');

epi = 1;
czx0 = ((1-lx)/(bx0(1)+bx1(1)))*tf([1 -ax1(1)],[epi -lx*epi 0 0 0 0 0 0 -(1-lx)],Ts,'variable','z^-1');

%Vertical SISO Controller
qzy = minreal(tf([1-ly],[1 -ly],Ts,'variable','z^-1')*tf([1 -ay1(1)],[by0(1)+by1(1)],Ts,'variable','z^-1'));
epi = 1/(1-ep);
czy = ((1-ly)/(by0(1)+by1(1)))*tf([1 -ay1(1)],[epi -ly*epi 0 0 0 0 0 0 -(1-ly)],Ts,'variable','z^-1');
epi = 1;
czy0 = ((1-ly)/(by0(1)+by1(1)))*tf([1 -ay1(1)],[epi -ly*epi 0 0 0 0 0 0 -(1-ly)],Ts,'variable','z^-1');

%%
si = 0.1;
ki = si^2/(1+si^2);

lzx = ki*czx*gzx(1);
lzy = ki*czy*gzy(1);

lzx0 = ki*czx0*gzx(1);
lzy0 = ki*czy0*gzy(1);

figure;
nyquist(lzx0); hold on;
nyquist(lzx, 'r--'); legend({'$\epsilon=0$','$\epsilon=10^{-3}$'},'Interpreter','Latex')
xlim([-1.5,1.5]); ylim([-1.5,1.5]); grid on; axis equal;

%%
si = 0.0299;
kc = si/(1+si^2);
mx = kc*czx*gzx(1)/(1+si*kc*czx*gzx(1));
mx0 = kc*czx0*gzx(1)/(1+si*kc*czx0*gzx(1));

bode_opt = bodeoptions();
bode_opt.FreqUnits = 'Hz';
bode_opt.XLim = [0.01, 100];
figure;
bodemag(1/mx0,'r-',1/mx,'b--',bode_opt); grid on;
legend({'$\epsilon=0$','$\epsilon=10^{-3}$'},'Interpreter','Latex')
%ylim([-30 0])

%% Full system

load(fname_RM);
RMorigx = Rmat(1).Data(:,:); % * 1e6; % RM(1) eq to RM(1,1)
[ny_x,nu_x] = size(RMorigx);
RMorigy = Rmat(4).Data(:,:); % * 1e6; % RM(4) eq to RM(2,2)
[ny_y,nu_y] = size(RMorigy);
assert(ny_x == ny_y);
assert(nu_x == nu_y);
[TOT_BPM,TOT_CM] = size(RMorigx);


R = RMorigx;
[Ux,Sx,Vx] = svd(R, 'econ');
MU =1.0*eye(size(R,2));
Dx = Sx.^2 / (Sx.^2+MU);
Ex = Sx / (Sx.^2+MU);
ISx = diag(diag(1./Sx));
IRMpx = Vx*Ex*Ux';

w_Hz = logspace(-2,2,400);
nw = length(w_Hz);
cx_w = freqresp(czx, w_Hz, 'Hz'); cx_w = cx_w(:);
cx0_w = freqresp(czx0, w_Hz, 'Hz'); cx0_w = cx0_w(:);
gx_w = freqresp(gzx(1), w_Hz, 'Hz'); gx_w = gx_w(:);

lx_w = freqresp(czx*gzx(1), w_Hz, 'Hz'); lx_w = lx_w(:);
lx0_w = freqresp(czx0*gzx(1), w_Hz, 'Hz'); lx0_w = lx0_w(:);

CR = IRMpx*R;
I = size(R,2);
norm_M = zeros(nw,1);
norm_M0 = zeros(nw,1);
for i=1:nw
    le = lx_w(i);
    l0 = lx0_w(i);
    M = diag(diag(Ex)*le)./(1+diag(Ex).*diag(Sx)*le);
    M0 = diag(diag(Ex)*l0)./(1+diag(Ex).*diag(Sx)*l0);
    
    %M = (I+CR*cx_w(i)*gx_w(i))\IRMpx*cx_w(i)*gx_w(i);
    %M0 = (I+CR*cx0_w(i)*gx_w(i))\IRMpx*cx0_w(i)*gx_w(i);
    norm_M(i) = norm(M,2);
    norm_M0(i) = norm(M0,2);
end

figure;
semilogx(w_Hz,20*log10(1./norm_M0),'r-',w_Hz,20*log10(1./norm_M),'b--'); grid on;
legend({'$\epsilon=0$','$\epsilon=10^{-3}$'},'Interpreter','Latex')
xlabel('Frequency [Hz]','Interpreter','Latex'); ylabel('Magnitude [dB]','Interpreter','Latex');
title('Max. $\vert\vert\Delta(\omega)\vert\vert_2$','Interpreter','Latex')

figure;
semilogx(w_Hz,20*log10(1./norm_M0)); grid on;




