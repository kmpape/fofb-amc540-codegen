function [Ao_x, Bo_x, Co_x, Ap_x, Bp_x, Cp_x, Ad_x, Cd_x,...
          Kfd_x, Kfx_x, Kcx_x, Kcd_x, P_x, Rlqr_x, Qlqr_x,...
          Ao_y, Bo_y, Co_y, Ap_y, Bp_y, Cp_y, Ad_y, Cd_y,...
          Kfd_y, Kfx_y, Kcx_y, Kcd_y, P_y, Rlqr_y, Qlqr_y] =...
          observer_regulator(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,nd,fname,print_msg)
    if nargin < 8; fname = ''; end
    if nargin < 9; print_msg = false; end
    %% Configure Diamond-I Storage Ring
    [~,TOT_CM_X] = size(RMorigx);
    [~,TOT_CM_Y] = size(RMorigy);
    ny_x = length(id_to_bpm_x);
    ny_y = length(id_to_bpm_y);
    nu_x = length(id_to_cm_x);
    nu_y = length(id_to_cm_y);

    %% Diamond-I Actuator Model
    aIx = 2*pi*500; % [rad/s]
    aIy = 2*pi*700; % [rad/s]
    Fs = 10*10^3; % sample frequency [Hz]
    Ts = 1/Fs; % sample time [s]
    %nd = 9; % number of delay samples [-]
    
    %% Frequencies for Sensitivity
    w_Hz = logspace(-2,4,300);
    w_Hz(w_Hz >= 5000) = [];
    
    %% Produce Gains for X & Y directions
    mu = 1;
    z = tf('z', Ts);
    dirs = {'X','Y'};
    for i=1:2
        if i == 1
            R = RMorigx(id_to_bpm_x, id_to_cm_x);
            aI = aIx;
        else
            R = RMorigy(id_to_bpm_y, id_to_cm_y);
            aI = aIy;
        end
        % Gains are generated in mode space
        [ny, ~] = size(R);
        [UR,SR,VR] = svd(R,'econ');

        bw = 1/(nd*Ts); 
        b = exp(-bw*Ts);
        a = exp(-aI*Ts);

        Kxtilde = zeros(ny,1); % x refers to state x here
        Kdtilde = zeros(ny,1);
        Rxtilde = zeros(ny,1);
        Ptilde = zeros(ny,1);
        sensitivity_lqr = zeros(length(w_Hz),ny);
        sensitivity_imc = zeros(length(w_Hz),ny);
        poles_imc = zeros(ny,1);
        for imode = 1:ny
            I = eye(nd+2);
            sval = SR(imode,imode);
            g_mp = sval*(1-a)/(z-a);
            g = g_mp * z^(-nd);
            c_imc = (1-b)/(1-a)*(z-a)*z^nd/(z^(nd+1)-b*z^nd-(1-b))*sval/(sval^2+mu);
            Smodei = inv(1+minreal(g*c_imc));
            pi_imc = pole(Smodei); % need to select the right pole
            pi_imc = pi_imc((imag(pi_imc)==0) & (real(pi_imc)>0));
            assert(length(pi_imc) == 1);
            if print_msg
                fprintf("%s, IMC closed loop sensitivity pole = %.4f Hz (mode = %d)\n",...
                    dirs{i}, -log(pi_imc)/2/pi/Ts, imode);
            end
            poles_imc(imode) = -log(pi_imc)/2/pi/Ts;
            Kdtilde(imode) = 1-pi_imc;

            % State-Space with delay
            Ap = [a,zeros(1,nd); eye(nd), zeros(nd,1)]; 
            Bp = [1-a; zeros(nd,1)]; 
            % Cp = [zeros(1,nd), sval]; % we don't use this. same weights
            % for every mode

            % [K,S,E] = dlqr(A,B,Q,R,N)
            % u[n] = -Kx[n]
            % J = Sum {x'Qx + u'Ru + 2*x'Nu}
            % x[n+1] = Ax[n] + Bu[n]
            Qlqr = diag([1,zeros(1,nd)]);
            Rlqr = sqrt((mu+sval^2)/sval^2);
            if false
                [Kc,Plqr,~] = dlqr(Ap,Bp,Qlqr,Rlqr,0);
                Kxtilde(imode) = Kc(1);
            else
                [Kc,Plqr,~] = dlqr(a,1-a,1,Rlqr,0);
                Kxtilde(imode) = Kc;
            end
            Rxtilde(imode) = Rlqr;
            Ptilde(imode) = Plqr(1,1);
            
            Ao = blkdiag(Ap,1);
            Bo = [Bp; 0];
            Co = [zeros(1,nd),sval, 1];
            Kf = [zeros(nd+1,1); Kdtilde(imode)];
            Kc = [Kxtilde(imode),zeros(1,nd),(1+Kxtilde(imode))/sval];
            I = eye(nd+2);
            c_lqr = ss((Ao-Bo*Kc)*(I-Kf*Co),(Ao-Bo*Kc)*Kf,Kc*(I-Kf*Co),Kc*Kf,Ts);
            OL = series(g,c_lqr);
            CL = 1/(1+OL);
            CL_w = freqresp(CL,w_Hz,'Hz'); CL_w = abs(CL_w(:));
            sensitivity_lqr(:,imode) = CL_w;
            tmp = freqresp(1/(1+g*c_imc),w_Hz,'Hz'); tmp = abs(tmp(:));
            sensitivity_imc(:,imode) = tmp;
            
            if false
                figure; bode(OL); grid on;
                figure; 
                semilogx(w_Hz,20*log10(CL_w),'b'); hold on;
                semilogx(w_Hz,20*log10(tmp),'r');
                grid on;
            end
        end

        if i==1
            Kfx_x = zeros(nu_x,ny_x);
            Kfd_x = UR*diag(Kdtilde)*UR';
            Kcx_x = VR*diag(Kxtilde)*VR';
            Kcd_x = VR*(SR\(eye(ny)+diag(Kxtilde)))*UR';
            P_x = VR*diag(Ptilde)*VR';
            Rlqr_x = VR*diag(Rxtilde)*VR';
            Qlqr_x = eye(nu_x);
            Slqr_x = sensitivity_lqr;
            Simc_x = sensitivity_imc;
            poles_imc_x = poles_imc;
        else
            Kfx_y = zeros(nu_y,ny_y);
            Kfd_y = UR*diag(Kdtilde)*UR';
            Kcx_y = VR*diag(Kxtilde)*VR';
            Kcd_y = VR*(SR\(eye(ny)+diag(Kxtilde)))*UR';
            P_y = VR*diag(Ptilde)*VR';
            Rlqr_y = VR*diag(Rxtilde)*VR';
            Qlqr_y = eye(nu_y);
            Slqr_y = sensitivity_lqr;
            Simc_y = sensitivity_imc;
            poles_imc_y = poles_imc;
        end
    end

    % Plants for observer and simulations (without delay)
    a_x = exp(-aIx*Ts);
    Ao_x = a_x*eye(nu_x);
    Bo_x = (1-a_x)*eye(nu_x);
    Co_x = RMorigx(id_to_bpm_x,id_to_cm_x);
    Ap_x = a_x*eye(TOT_CM_X);
    Bp_x = (1-a_x)*eye(TOT_CM_X);
    Cp_x = RMorigx;

    a_y = exp(-aIy*Ts);
    Ao_y = a_y*eye(nu_y);
    Bo_y = (1-a_y)*eye(nu_y);
    Co_y = RMorigy(id_to_bpm_y,id_to_cm_y);
    Ap_y = a_y*eye(TOT_CM_Y);
    Bp_y = (1-a_y)*eye(TOT_CM_Y);
    Cp_y = RMorigy;

    % Disturbance
    Ad_x = eye(ny_x);
    Cd_x = eye(ny_x);
    Ad_y = eye(ny_y);
    Cd_y = eye(ny_y);

    if ~isempty(fname)
       save(fname, 'Ao_x', 'Bo_x', 'Co_x', 'Ap_x', 'Bp_x', 'Cp_x', 'Ad_x', 'Cd_x',...
                   'Kfd_x', 'Kfx_x', 'Kcx_x', 'Kcd_x', 'P_x', 'Rlqr_x', 'Qlqr_x',...
                   'Ao_y', 'Bo_y', 'Co_y', 'Ap_y', 'Bp_y', 'Cp_y', 'Ad_y', 'Cd_y',...
                   'Kfd_y', 'Kfx_y', 'Kcx_y', 'Kcd_y', 'P_y', 'Rlqr_y', 'Qlqr_y');
       save(strrep(fname,'.mat','_analysis.mat'), 'poles_imc_x', 'poles_imc_y',...
            'w_Hz', 'Slqr_x', 'Simc_x', 'Slqr_y', 'Simc_y');
    end
end



