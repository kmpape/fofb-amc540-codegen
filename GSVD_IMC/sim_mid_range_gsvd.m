function [y_sim, u_sim, us_sim, uf_sim] = sim_mid_range_gsvd(...
            n_delay, doff,...
            A, B, C, D,...                   % Plant
            Ks, Acs, Bcs, Ccs, Dcs,...       % Slow controller
            Kf, Acf, Bcf, Ccf, Dcf,...       % Fast controller
            Ps, Ams, Bms, Cms, Dms,...       % Slow plant model
            Pf, Amf, Bmf, Cmf, Dmf,...       % Fast plant model
            id_to_bpm, slow_to_id, fast_to_id, ol_mode, network_scaling)
n_samples = size(doff, 2);
% assert(norm(Dcs) == 0);
% assert(norm(Dcf) == 0);
% assert(norm(Dms) == 0);
% assert(norm(Dmf) == 0);
use_single = true;

x_new   = double(zeros(size(A,1), 1));   x_old = x_new;
xcs_new = double(zeros(size(Acs,1), 1)); xcs_old = xcs_new;
xcf_new = double(zeros(size(Acf,1), 1)); xcf_old = xcf_new;
xms_new = double(zeros(size(Ams,1), 1)); xms_old = xms_new;
yms = double(zeros(size(Cms,1),1));
xmf_new = double(zeros(size(Amf,1), 1)); xmf_old = xmf_new;
ymf = double(zeros(size(Cms,1),1));

y_sim = double(zeros(size(C,1), n_samples));
u_sim = double(zeros(size(B,2), n_samples));
us_sim = double(zeros(size(Ccs,1), n_samples));
uf_sim = double(zeros(size(Ccf,1), n_samples));

Acs = double(Acs); Bcs = double(Bcs); Ccs = double(Ccs); Dcs = double(Dcs);
Acf = double(Acf); Bcf = double(Bcf); Ccf = double(Ccf); Dcf = double(Dcf);
Ams = double(Ams); Bms = double(Bms); Cms = double(Cms); Dms = double(Dms);
Amf = double(Amf); Bmf = double(Bmf); Cmf = double(Cmf); Dmf = double(Dmf);

if use_single == true
    Ks = single(Ks);
    Kf = single(Kf);
    Ps = single(Ps);
    Pf = single(Pf);
else
    Ks = double(Ks);
    Kf = double(Kf);
    Ps = double(Ps);
    Pf = double(Pf);
end

truncate_to_17_bits = false;
if truncate_to_17_bits
    tmp = max(max(abs(Ks)));
    Ks = double(round(Ks/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Acs)));
    Acs = double(round(Acs/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Bcs)));
    Bcs = double(round(Bcs/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Ccs)));
    Ccs = double(round(Ccs/tmp*2^17,0)/2^17)*tmp; 
    Dcs = double(round(Dcs*2^17,0)/2^17);
    tmp = max(max(abs(Kf)));
    Kf = double(round(Kf/tmp*2^17,0)/2^17);
    tmp = max(max(abs(Acf)));
    Acf = double(round(Acf/tmp*2^17,0)/2^17)*tmp;
    tmp = max(max(abs(Bcf)));
    Bcf = double(round(Bcf/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Ccf)));
    Ccf = double(round(Ccf/tmp*2^17,0)/2^17)*tmp; 
    Dcf = double(round(Dcf*2^17,0)/2^17);
    tmp = max(max(abs(Ps)));
    Ps = double(round(Ps/tmp*2^17,0)/2^17)*tmp;
    tmp = max(max(abs(Ams)));    
    Ams = double(round(Ams/tmp*2^17,0)/2^17)*tmp;
    tmp = max(max(abs(Bms)));
    Bms = double(round(Bms/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Cms)));
    Cms = double(round(Cms/tmp*2^17,0)/2^17)*tmp; 
    Dms = double(round(Dms*2^17,0)/2^17);
    tmp = max(max(abs(Pf)));
    Pf = double(round(Pf/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Amf)));
    Amf = double(round(Amf/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Bmf)));
    Bmf = double(round(Bmf/tmp*2^17,0)/2^17)*tmp; 
    tmp = max(max(abs(Cmf)));
    Cmf = double(round(Cmf/tmp*2^17,0)/2^17)*tmp; 
    Dmf = double(round(Dmf*2^17,0)/2^17);
end

%%
a = 1e3/1e6;

% XCS_MAX = double(5 / (abs(Ccs(1,1))) * ones(size(Acs,1), 1));
% XCF_MAX = double(5 / (abs(Ccf(1,1))) * ones(size(Acf,1), 1));

do_print = false;
for ii = 1 : n_samples
    if ol_mode
        y_sim(:, ii) = doff(:, ii);
    else
        y_sim(:, ii) = C * x_old + D * double(u_sim(:, ii)) / a + doff(:, ii);
    end
    
    if ii > n_delay
        
        if use_single == true
            fb_signal = single(y_sim(id_to_bpm, ii-n_delay)*1000 + yms + ymf);
        else
            fb_signal = y_sim(id_to_bpm, ii-n_delay)*1000 + yms + ymf;
        end

        % control inputs
        in_xcs = double(Ks * fb_signal);
        if do_print; fprintf("in_xcs=[%.12f,%.12f,%.12f,]\n",in_xcs(1),in_xcs(2),in_xcs(3)); end
        xcs_new = Acs * xcs_old + Bcs * in_xcs;
        us_sim(:,ii) = Ccs * xcs_old + Dcs * in_xcs;
        if do_print; fprintf("us_sim=[%.12f,%.12f,%.12f,]\n",us_sim(1,ii),us_sim(2,ii),us_sim(3,ii)); end
        xcs_old = xcs_new;

        in_xcf = double(Kf * fb_signal);
        if do_print; fprintf("in_xcf=[%.12f,%.12f,%.12f,]\n",in_xcf(1),in_xcf(2),in_xcf(3)); end
        xcf_new = Acf * xcf_old + Bcf * in_xcf;
        uf_sim(:,ii) = Ccf * xcf_old + Dcf * in_xcf;
        % xcf_new = max(-XCF_MAX, min(Acf * xcf_old + Bcf * in_xcf, XCF_MAX));
        % uf_sim(:,ii) = Ccf * xcf_old + Dcf * in_xcf;
        if do_print; fprintf("uf_sim=[%.12f,%.12f,%.12f,]\n",uf_sim(1,ii),uf_sim(2,ii),uf_sim(3,ii)); end
        xcf_old = xcf_new;

        % model outputs
        if use_single == true
            in_xms = double(Ps * single(us_sim(:,ii)));
        else
            in_xms = double(Ps * us_sim(:,ii));
        end
        if do_print; fprintf("in_xms=[%.12f,%.12f,%.12f,]\n",in_xms(1),in_xms(2),in_xms(3)); end
        xms_new = Ams * xms_old + Bms * in_xms;
        yms = Cms * xms_old + Dms * in_xms;
        if do_print; fprintf("yms=[%.12f,%.12f,%.12f,]\n",yms(1),yms(2),yms(3)); end
        xms_old = xms_new;

        if use_single == true
            in_xmf = double(Pf * single(uf_sim(:,ii)));
        else
            in_xmf = double(Pf * uf_sim(:,ii));
        end
        if do_print; fprintf("in_xmf=[%.12f,%.12f,%.12f,]\n",in_xmf(1),in_xmf(2),in_xmf(3)); end
        xmf_new = Amf * xmf_old + Bmf * in_xmf;
        ymf = Cmf * xmf_old + Dmf * in_xmf;
        if do_print; fprintf("ymf=[%.12f,%.12f,%.12f,]\n",ymf(1),ymf(2),ymf(3)); end
        xmf_old = xmf_new;

        u_sim(slow_to_id, ii+1) = us_sim(:,ii);
        u_sim(fast_to_id, ii+1) = uf_sim(:,ii);
        
        u_sim(:, ii+1)= round(network_scaling*u_sim(:, ii+1),0)/network_scaling;
    end
    
    if ~ol_mode
        x_new = A * x_old + B * double(u_sim(:, ii+1)) / a;
        x_old = x_new;
    end
end

y_sim = y_sim';
u_sim = u_sim';
us_sim = us_sim';
uf_sim = uf_sim';
end
