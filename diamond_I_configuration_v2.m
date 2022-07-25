function [id_to_bpm_x,  id_to_bpm_y,...
          slow_to_id_x, slow_to_id_y,...
          fast_to_id_x, fast_to_id_y] = diamond_I_configuration_v2(RMorigx,RMorigy)
%%
% RMorigx = Rmat(1).Data(:,:);
% RMorigy = Rmat(4).Data(:,:);
[TOT_BPM,TOT_CM] = size(RMorigx);

bad_bpm = sort([76, 79],'asc');
bad_cm = sort([],'asc');
all_bpm = 1 : 1 : TOT_BPM;
all_cm = 1 : 1 : TOT_CM;
cell_cm_i =  [1,8,16,23,30,37,44,51,58,67,74,80,87,96,103,110,117,124,131,138,145,152,159,166,TOT_CM];
cell_bpm_i = [1,8,16,23,30,37,44,51,58,67,74,81,88,97,104,111,118,125,132,139,146,153,160,167,TOT_BPM];
assert(length(cell_cm_i)-1==24);
assert(length(cell_bpm_i)-1==24);

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
nu_good = size(RMx,2);

num_slow = 4*ones(1,24);
num_slow(1:2) = 5;
num_slow(9:10) = 5;
num_slow(12:14) = 5;
num_slow(23:24) = 5;

num_fast = 2*ones(1,24);
num_fast(9) = 3;
num_fast(13) = 3;
num_fast(3) = 3;
num_fast(22) = 3;
num_fast(21) = 3;


slow_inds_x = []; fast_inds_x = [];
slow_inds_y = []; fast_inds_y = [];
orm_vecnorm_x = vecnorm(RMx,2,1);
orm_vecnorm_y = vecnorm(RMy,2,1);
for i=1:24
    if (i<24); inds = cell_cm_i(i):1:cell_cm_i(i+1)-1; else; inds = cell_cm_i(i):1:cell_cm_i(end); end
    if length(inds) < num_slow(i)+num_fast(i)
        i
        break
    end
    notinds = setdiff(used_cm,inds);
    tmp = orm_vecnorm_x;
    tmp(notinds) = 0;
    [~,indssorted] = sort(tmp,'desc');
    slow_inds_x = [slow_inds_x, indssorted(1:num_slow(i))];
    fast_inds_x = [fast_inds_x, indssorted(num_slow(i)+1:num_slow(i)+num_fast(i))];
    tmp = orm_vecnorm_y;
    tmp(notinds) = 0;
    [~,indssorted] = sort(tmp,'desc');
    slow_inds_y = [slow_inds_y, indssorted(1:num_slow(i))];
    fast_inds_y = [fast_inds_y, indssorted(num_slow(i)+1:num_slow(i)+num_fast(i))];
end
slow_inds_x = sort(slow_inds_x,'asc');
fast_inds_x = sort(fast_inds_x,'asc');
slow_inds_y = sort(slow_inds_y,'asc');
fast_inds_y = sort(fast_inds_y,'asc');

if 0
    figure;
    subplot(2,1,1);
    plot(orm_vecnorm_x,'k'); hold on;
    plot(slow_inds_x,orm_vecnorm_x(slow_inds_x),'b','Linestyle','None','Marker','x');
    plot(fast_inds_x,orm_vecnorm_x(fast_inds_x),'r','Linestyle','None','Marker','o');
    title('X CM Selection (x=slow,o=fast)')
    xlim([1 TOT_CM])
    for i=1:24; plot([cell_cm_i(i),cell_cm_i(i)],[0 30],'k--'); end
    subplot(2,1,2);
    plot(orm_vecnorm_y,'k'); hold on;
    plot(slow_inds_y,orm_vecnorm_y(slow_inds_y),'b','Linestyle','None','Marker','x');
    plot(fast_inds_y,orm_vecnorm_y(fast_inds_y),'r','Linestyle','None','Marker','o');
    for i=1:24; plot([cell_cm_i(i),cell_cm_i(i)],[0 30],'k--'); end
    title('Y CM Selection (x=slow,o=fast)')
    xlim([1 TOT_CM])
end

is_slow_x = zeros(nu_good,1); is_slow_x(slow_inds_x) = 1;
is_fast_x = zeros(nu_good,1); is_fast_x(fast_inds_x) = 1;
is_slow_y = zeros(nu_good,1); is_slow_y(slow_inds_y) = 1;
is_fast_y = zeros(nu_good,1); is_fast_y(fast_inds_y) = 1;

slow_to_id_x = slow_inds_x;
fast_to_id_x = fast_inds_x;
slow_to_id_y = slow_inds_y;
fast_to_id_y = fast_inds_y;

assert(sum((is_slow_x == is_fast_x) & (is_slow_x > 0))==0);
assert(sum((is_slow_y == is_fast_y) & (is_slow_y > 0))==0);
assert(isempty(bad_cm));

ns = sum(is_slow_x);
nf = sum(is_fast_x);
nu = ns + nf;
ny = ns;

[~, active_bpms_x] = licols(RMx(:,is_slow_x==1)');
Rx = RMx(active_bpms_x,:);
Rsx = Rx(:,is_slow_x==1);
Rfx = Rx(:,is_fast_x==1);
assert(rank(Rsx)==ns);

[~, active_bpms_y] = licols(RMy(:,is_slow_y==1)');
Ry = RMy(active_bpms_y,:);
Rsy = Ry(:,is_slow_y==1);
Rfy = Ry(:,is_fast_y==1);
assert(rank(Rsy)==ns);

id_to_bpm_x = used_bpm(active_bpms_x);
id_to_bpm_y = used_bpm(active_bpms_y);

is_bpm_x = zeros(1,length(used_bpm)); is_bpm_x(active_bpms_x) = 1;
is_bpm_y = zeros(1,length(used_bpm)); is_bpm_y(active_bpms_y) = 1;


if 0
    figure;
    subplot(2,1,1);
    plot(is_bpm_x); hold on;
    plot(active_bpms_x,is_bpm_x(active_bpms_x),'r','Linestyle','None','Marker','x');
    title('X BPM Selection'); xlim([1 length(used_bpm)]);
    for i=1:24; plot([cell_bpm_i(i),cell_bpm_i(i)],[0 1],'k--'); end
    subplot(2,1,2);
    plot(is_bpm_y); hold on;
    plot(active_bpms_y,is_bpm_y(active_bpms_y),'r','Linestyle','None','Marker','x');
    title('X BPM Selection'); xlim([1 length(used_bpm)]);
    for i=1:24; plot([cell_bpm_i(i),cell_bpm_i(i)],[0 1],'k--'); end
    title('Y BPM Selection')
    xlim([1 length(used_bpm)])
end


end


