function [id_to_bpm_x, id_to_cm_x, id_to_bpm_y, id_to_cm_y] = diamond_I_configuration_v4(RMorigx,RMorigy)
%%
TOT_BPM = 173;
TOT_CM = 172;

% This is the configuration implemented on the 03.08.2022 at DLS
bad_bpm_x = [32  75  78 120 155]+1;
bad_bpm_y = [18  25  62  75  78  83 141]+1;
bad_cm_x = [76]+1;
bad_cm_y = [76]+1;

% Generate IDs as used by the FOFB (note: codegen will subtract 1)
id_to_bpm_x = 1:1:TOT_BPM;
id_to_bpm_x(bad_bpm_x) = [];
id_to_bpm_y = 1:1:TOT_BPM;
id_to_bpm_y(bad_bpm_y) = [];
id_to_cm_x = 1:1:TOT_CM;
id_to_cm_x(bad_cm_x) = [];
id_to_cm_y = 1:1:TOT_CM;
id_to_cm_y(bad_cm_y) = [];

% Make sure we still have full rank
assert(rank(RMorigx(id_to_bpm_x, id_to_cm_x)) == length(id_to_bpm_x));
assert(rank(RMorigy(id_to_bpm_y, id_to_cm_y)) == length(id_to_bpm_y));

for i = 1:length(bad_bpm_x); assert(any(bad_bpm_x(i)==id_to_bpm_x)==0); end
for i = 1:length(bad_cm_x); assert(any(bad_cm_x(i)==id_to_cm_x)==0); end
for i = 1:length(bad_bpm_y); assert(any(bad_bpm_y(i)==id_to_bpm_y)==0); end
for i = 1:length(bad_cm_y); assert(any(bad_cm_y(i)==id_to_cm_y)==0); end

end


