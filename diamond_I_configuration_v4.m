function [id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v4(RMorigx,RMorigy,square_config)
if nargin < 3
    square_config = false;
end

%%
TOT_BPM = 173;
TOT_CM = 172;

% This is the configuration implemented on the 03.08.2022 at DLS
if ~square_config
    bad_bpm_x = [32  75  78 120 155]+1;
    bad_bpm_y = [18  25  62  75  78  83 141]+1;
    bad_cm_x = [76]+1;
    bad_cm_y = [76]+1;
else
    bad_bpm_x = [32  75  78 120 155]+1;
    bad_bpm_y = [18  25  62  75  78  83 141]+1;
    bad_cm_x = sort([76,58,155,120],'asc')+1;
    bad_cm_y = sort([76,58,141,18,25,83],'asc')+1;    
    bad_cm_x = [76 32 120 155]+1;
    bad_cm_y = [76 18 25 62 83 141]+1;
end

% Generate IDs as used by the FOFB (note: codegen will subtract 1)
id_to_bpm_x = 1:1:TOT_BPM;
id_to_bpm_x(bad_bpm_x) = [];
id_to_bpm_y = 1:1:TOT_BPM;
id_to_bpm_y(bad_bpm_y) = [];
id_to_cm_x = 1:1:TOT_CM;
id_to_cm_x(bad_cm_x) = [];
id_to_cm_y = 1:1:TOT_CM;
id_to_cm_y(bad_cm_y) = [];

% fprintf("BPM x=%d, CM x=%d, BPM y=%d, CM y=%d\n",...
%     length(id_to_bpm_x),length(id_to_cm_x),length(id_to_bpm_y),length(id_to_cm_y))

% Make sure we still have full rank
assert(rank(RMorigx(id_to_bpm_x, id_to_cm_x)) == length(id_to_bpm_x));
assert(rank(RMorigy(id_to_bpm_y, id_to_cm_y)) == length(id_to_bpm_y));

for i = 1:length(bad_bpm_x); assert(any(bad_bpm_x(i)==id_to_bpm_x)==0); end
for i = 1:length(bad_cm_x); assert(any(bad_cm_x(i)==id_to_cm_x)==0); end
for i = 1:length(bad_bpm_y); assert(any(bad_bpm_y(i)==id_to_bpm_y)==0); end
for i = 1:length(bad_cm_y); assert(any(bad_cm_y(i)==id_to_cm_y)==0); end

%%
end


