function [y_sim, u_sim] = sim_standard_imc(...
    nsamples, n_delay, id_to_bpm, id_to_cm, dist,...
    Ac, Bc, Cc, Dc,... % CONTROLLER STATE-SPACE
    A, B, C, D, network_scaling, open_loop, run_ntimes) % PLANT STATE-SPACE


% Input to controller:      NANOMETERS
% Output of controller:     AMPERES
% DSP scaling out:          1e6
% VME scaling in:           -1/1e6 (fixed)
if ~exist('open_loop', 'var') open_loop = false; end
if ~exist('run_ntimes', 'var') run_ntimes = 1; end

%%
use_single = true;

% Does not support feedthrough on plant
assert(sum(sum(abs(D))) == 0);
assert(size(dist, 2) == nsamples);
%assert(size(dist, 1) == size(C, 1));

% Dimensions Plant == Dimension Model
[nx, ~] = size(A);
[ny, nu] = size(D);

% Dimensions Controller
[nxc, ~] = size(Ac);
[nyc, ~] = size(Cc);
%assert(nyc == nu);
if use_single == true
    x_ctr_old = single(zeros(nxc,1));
    Ac = single(Ac);
    Bc = single(Bc);
    Cc = single(Cc);
    Dc = single(Dc);
else
    x_ctr_old = zeros(nxc,1);
end

% Variables Plant
x_sim_new = zeros(nx, 1);
x_sim_old = zeros(nx, 1);
y_sim = zeros(ny, nsamples);

% Variables Controller
u_sim = zeros(nu, nsamples);

%%
for mm = 1 : run_ntimes
    for ii = 1 : nsamples

        % Measure Plant Output
        if open_loop == false
            y_sim(:, ii) = C * x_sim_new + dist(:, ii);
        else
            y_sim(:, ii) = dist(:, ii);
        end

        % Update Controller
        if (ii > n_delay) || (mm > 1)
            ind_ii = ii-n_delay;
            if ind_ii <= 0
                ind_ii = nsamples+ind_ii;
            end
            if use_single == true
                y_meas = single(round(y_sim(id_to_bpm, ind_ii)*1000,0));
            else
                y_meas = double(round(y_sim(id_to_bpm, ind_ii)*1000,0));
            end

            x_ctr_new = Ac*x_ctr_old + Bc*y_meas;
            u_ctr = Cc*x_ctr_old + Dc*y_meas;
            x_ctr_old = x_ctr_new;

            u_sim(id_to_cm, ii) = double(round(u_ctr*network_scaling,0)/(network_scaling));
        end

        % Update Plant Model
        x_sim_new = A * x_sim_old + B * u_sim(:, ii)*1000;
        x_sim_old = x_sim_new;
    end
end

%%
y_sim = y_sim';
u_sim = u_sim';


end
