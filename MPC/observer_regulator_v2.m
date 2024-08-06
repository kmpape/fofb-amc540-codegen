function [Aoi_x, Boi_x, Coi_x, Kfi_x, Kci_x, ...
          Aoi_y, Boi_y, Coi_y, Kfi_y, Kci_y] =...
          observer_regulator_v2(RMorigx,RMorigy,id_to_bpm_x, id_to_cm_x, id_to_bpm_y,id_to_cm_y,nd,fname,print_msg)
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
        Aoi = zeros(nd+2,nd+2,ny);
        Boi = zeros(nd+2,1,ny);
        Coi = zeros(1,nd+2,ny);
        Kfi = zeros(nd+2,1,ny);
        Kci = zeros(1,nd+2,ny);
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
            tmp = 0.8*poles_imc(imode);
            %Kdtilde(imode) = 1-pi_imc;
            Kdtilde(imode) = 1-exp(-tmp*2*pi*Ts);

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
            Aoi(:,:,imode) = Ao;
            Boi(:,:,imode) = Bo;
            Coi(:,:,imode) = Co;
            Kfi(:,:,imode) = Kf;
            Kci(:,:,imode) = Kc;
        end

        if i==1
            Aoi_x = Aoi;
            Boi_x = Boi;
            Coi_x = Coi;
            Kfi_x = Kfi;
            Kci_x = Kci;
        else
            Aoi_y = Aoi;
            Boi_y = Boi;
            Coi_y = Coi;
            Kfi_y = Kfi;
            Kci_y = Kci;
        end
    end
    if ~isempty(fname)
       save(fname, 'Aoi_x', 'Boi_x', 'Coi_x', 'Kfi_x', 'Kci_x', ...
                   'Aoi_y', 'Boi_y', 'Coi_y', 'Kfi_y', 'Kci_y');
    end
end



