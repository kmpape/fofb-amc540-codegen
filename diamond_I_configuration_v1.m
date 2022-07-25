function [id_to_bpm, slow_to_id, fast_to_id] = diamond_I_configuration_v1(RMorigx,RMorigy)
%%
% RMorigx = Rmat(1).Data(:,:);
% RMorigy = Rmat(4).Data(:,:);
[TOT_BPM,TOT_CM] = size(RMorigx);

bad_bpm = sort([76, 79],'asc');
bad_cm = sort([],'asc');
cell_cm_i =  [1,8,16,23,30,37,44,51,58,67,74,80,87,96,103,110,117,124,131,138,145,152,159,166,TOT_CM];
cell_bpm_i = [1,8,16,23,30,37,44,51,58,67,74,81,88,97,104,111,118,125,132,139,146,153,160,167,TOT_BPM];
assert(length(cell_cm_i)-1==24);
assert(length(cell_bpm_i)-1==24);
assert(isempty(bad_cm));

% Remove unused BPM/CM
used_bpm = 1 : 1 : TOT_BPM;
used_cm = 1 : 1 : TOT_CM;
RMx = RMorigx; RMy = RMorigy;
RMx(:,bad_cm) = [];  RMy(:,bad_cm) = [];
RMx(bad_bpm,:) = []; RMy(bad_bpm,:) = [];
used_bpm(bad_bpm) = [];
used_cm(bad_cm) = [];
for i = 1 : length(bad_cm)
    cell_cm_i = cell_cm_i - (cell_cm_i>=bad_cm(i)-(i-1));
end
for i = 1 : length(bad_bpm)
    cell_bpm_i = cell_bpm_i - (cell_bpm_i>=bad_bpm(i)-(i-1));
end
nu_good = size(RMorigx,2);

% Diamond-I storage ring reconfiguration: ns=108=4.5*24 and nf=60=2.5*24
% Some cells have different number of magnets
special_cells = [2,24,9,11,13];
slow_pattern = containers.Map('keyType','char','ValueType','any');
slow_pattern('even') = [1,0,1,0,1,0,1]; slow_pattern('odd') = [1,0,1,1,1,0,1];
slow_pattern('2') = [1,1,1,0,0,0,0,1]; slow_pattern('24') = [1,0,0,0,1,1];
slow_pattern('9') = [1,0,1,1,1,1,0,0,1]; slow_pattern('11') = [1,0,1,1,0,1];
slow_pattern('13')= [1,0,1,0,1,1,0,1,1];

fast_pattern = containers.Map('keyType','char','ValueType','any');
fast_pattern('even') = [0,1,0,1,0,1,0]; fast_pattern('odd') = [0,1,0,0,0,1,0];
fast_pattern('2') = [0,0,0,1,1,1,0,0]; fast_pattern('24') = [0,1,1,1,0,0];
fast_pattern('9') = [0,1,0,0,0,0,1,0,0]; fast_pattern('11') = [0,1,0,0,1,0];
fast_pattern('13') = [0,1,0,0,0,0,1,0,0];

is_slow = zeros(nu_good,1);
is_fast = zeros(nu_good,1);
for i = 1 : 24
    if any(i == special_cells)
        is_slow(cell_cm_i(i):cell_cm_i(i+1)-1) = slow_pattern(char(string(i)));
        is_fast(cell_cm_i(i):cell_cm_i(i+1)-1) = fast_pattern(char(string(i)));
    elseif mod(i,2) == 0
        is_slow(cell_cm_i(i):cell_cm_i(i+1)-1) = slow_pattern('even');
        is_fast(cell_cm_i(i):cell_cm_i(i+1)-1) = fast_pattern('even');
    else
        is_slow(cell_cm_i(i):cell_cm_i(i+1)-1) = slow_pattern('odd');
        is_fast(cell_cm_i(i):cell_cm_i(i+1)-1) = fast_pattern('odd');
    end
end
assert(sum((is_slow == is_fast) & (is_slow > 0))==0);
ns = sum(is_slow);
nf = sum(is_fast);
ny = ns;

% Map from computed corrector output to network ID
slow_to_id = zeros(1, ns); k = 1;
for i = 1 : length(is_slow)
    if is_slow(i) == 1
        slow_to_id(k) = used_cm(i); k = k + 1;
    end
end
fast_to_id = zeros(1, nf); k = 1;
for i = 1 : length(is_fast)
    if is_fast(i) == 1
        fast_to_id(k) = used_cm(i); k = k + 1;
    end
end

% Need to select ns BPMs, so that rank(Rs)=ns and rank(Rf)=nf
[~, active_bpms] = licols(RMx(:,is_slow==1)');
Rx = RMx(active_bpms,:);
Rsx = Rx(:,is_slow==1);
assert(rank(Rsx)==ns);

Ry = RMy(active_bpms,:);
Rsy = Ry(:,is_slow==1);
assert(rank(Rsy)==ns);

% Map network id to BPM
id_to_bpm = zeros(1, ny);
for i = 1 : ny
    id_to_bpm(i) = used_bpm(active_bpms(i));
end


end


