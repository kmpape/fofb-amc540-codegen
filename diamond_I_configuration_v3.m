function [id_to_bpm, slow_to_id, fast_to_id] = diamond_I_configuration_v3(RMorigx,RMorigy)
%%
if (0)
    fname_RM = '/ORMS/GoldenBPMResp.mat';
    load(fname_RM);
    RMorigx = Rmat(1).Data(:,:);
    RMorigy = Rmat(4).Data(:,:);
end
[TOT_BPM,TOT_CM] = size(RMorigx);

bad_bpm = sort([76, 79],'asc');
bad_cm = sort([],'asc');
cell_cm_i =  [1,8,16,23,30,37,44,51,58,67,74,80,87,96,103,110,117,124,131,138,145,152,159,166,TOT_CM];
cell_bpm_i = [1,8,16,23,30,37,44,51,58,67,74,81,88,97,104,111,118,125,132,139,146,153,160,167,TOT_BPM];
assert(length(cell_cm_i)-1==24);
assert(length(cell_bpm_i)-1==24);
assert(isempty(bad_cm));

ns = 24*4;
nf = min(3*24,4*16); % would like to have 3*24 but > 4*16
ny = 24*4;

ebpm_pattern = containers.Map('keyType','char','ValueType','any');
slow_pattern = containers.Map('keyType','char','ValueType','any');
fast_pattern = containers.Map('keyType','char','ValueType','any');
ebpm_pattern('base') = [1,0,1,0,1,0,1];
ebpm_pattern('2')    = [1,0,1,0,0,1,0,1];
ebpm_pattern('9')    = [0,0,1,0,1,0,1,0,1];
ebpm_pattern('11')   = [1,1,0,0,1,0,1];
ebpm_pattern('13')   = [0,0,1,0,1,0,1,0,1];


slow_pattern('base') = [1,0,1,0,1,0,1];
slow_pattern('2')    = [1,0,1,0,0,1,0,1];
slow_pattern('9')    = [0,0,1,0,1,0,1,0,1];
slow_pattern('11')   = [1,0,1,1,0,1];
slow_pattern('13')   = [0,0,1,0,1,0,1,0,1];


fast_pattern('base') = [0,1,0,1,0,1,0];
fast_pattern('2')    = [0,1,0,0,1,0,1,0];
fast_pattern('9')    = [0,0,0,1,0,1,0,1,0];
fast_pattern('11')   = [0,1,0,0,1,0];
fast_pattern('13')   = [1,0,0,1,0,1,0,1,0];
% % need to have max 4*16 fast actuators
fast_pattern('1') = [0,1,0,0,0,1,0];
fast_pattern('3') = [0,1,0,0,0,1,0];
fast_pattern('5') = [0,1,0,0,0,1,0];
fast_pattern('7') = [0,1,0,0,0,1,0];
fast_pattern('15') = [0,1,0,0,0,1,0];
fast_pattern('17') = [0,1,0,0,0,1,0];
fast_pattern('19') = [0,1,0,0,0,1,0];
fast_pattern('21') = [0,1,0,0,0,1,0];

% ns = ns + 1;
% nf = nf - 1;
% ny = ny + 1;
% ebpm_pattern('16')   = [1,0,1,1,1,0,1]; %%%
% slow_pattern('16')   = [1,0,1,1,1,0,1]; %%%
% fast_pattern('16') = [0,1,0,0,0,1,0]; %%%

special_cells = [2,9,11,13];
special_cells_fast = [2,9,11,13,1,3,5,7,15,17,19,21];
is_ebpm = [];
is_slow = [];
is_fast = [];
for i = 1 : 24
    if any(i == special_cells)
        is_ebpm = [is_ebpm, ebpm_pattern(char(string(i)))];
        is_slow = [is_slow, slow_pattern(char(string(i)))];
    else
        is_ebpm = [is_ebpm, ebpm_pattern('base')];
        is_slow = [is_slow, slow_pattern('base')];
    end
    if any(i == special_cells_fast)
        is_fast = [is_fast, fast_pattern(char(string(i)))];
    else
        is_fast = [is_fast, fast_pattern('base')];
    end
end
assert(sum(is_slow) == ns);
assert(sum(is_fast) == nf);
assert(sum(is_ebpm) == ny);
assert(sum((is_slow==is_fast) & (is_slow==1)) == 0);
assert(rank(RMorigx(is_ebpm==1,is_slow==1)) == ns);
assert(rank(RMorigy(is_ebpm==1,is_slow==1)) == ns);

all_nums = 1:1:173;
id_to_bpm  = all_nums(is_ebpm==1);
all_nums = 1:1:172;
slow_to_id = all_nums(is_slow==1);
fast_to_id = all_nums(is_fast==1);

assert(rank(RMorigx(id_to_bpm,slow_to_id)) == ns);
assert(rank(RMorigy(id_to_bpm,slow_to_id)) == ns);

for i = 1:length(bad_bpm); assert(any(bad_bpm(i)==id_to_bpm)==0); end
for i = 1:length(bad_cm); assert(any(bad_cm(i)==slow_to_id)==0); end
for i = 1:length(bad_cm); assert(any(bad_cm(i)==fast_to_id)==0); end

% Code limitation
%assert(nf <= 4*16);

end


