function generate_filter_from_TF(foldername, vec_len, num_f, den_f, filter_name,...
    filter_base_type, print_test_data, datasection, alignment, integrator_max)
% 
% Generates C-code for a discrete-time IIR filter of the form
% y(z)/u(z) = (num_f(1)+num_f(2)*z^-1+...+num_f(N)*z^-N) /...
%             (den_f(1)+den_f(2)*z^-1+...+num_f(M)*z^-M)
% with M >= Z and where y and u are vectors of length vec_len.
%

% Example filter parameters
% vec_len = 23;
% num_f = [0.1, 0.2, 0.3];
% den_f = [1, -0.4, 0, -0.6];
% filter_name = '1st';
if nargin < 5
    filter_base_type = 'float';
end
if nargin < 6
print_test_data = 1;
end
if nargin < 7
    datasection = '';
end
if nargin < 8
    alignment = 16;
end

%% Preamble
assert(isscalar(vec_len));
assert(ischar(filter_name));
assert(ischar(filter_base_type));

% Remove trailing zeros in numerator and denumerator
den_f = den_f(1 : find(den_f, 1, 'last'));
num_f = num_f(1 : find(num_f, 1, 'last'));

% Check lengths
n_den = length(den_f);
n_num = length(num_f);
% assert(n_den >= n_num);

% Check base type
assert(strcmp(filter_base_type,'float') ||...
    strcmp(filter_base_type,'double'));

% normalize first output coefficient
if (abs(den_f(1)) > 0) && (den_f(1) ~= 1)
    num_f = num_f ./ den_f(1);
    den_f(2:end) = den_f(2:end) ./ den_f(1);
    den_f(1) = 1;
end

%foldername = sprintf('DTF_%s/', filter_name);
%status = mkdir(foldername);

filename_c = sprintf('DTF_%s.c', filter_name);
filename_h = sprintf('DTF_%s.h', filter_name);
prefix = sprintf('DTF_%s_', filter_name);
if print_test_data == 1
    filename_test_data = sprintf('DTF_%s_UNIT_TEST.h', filter_name);
end

% vector length name
flen = sprintf('%sLEN', prefix);
% array type name
ftype = sprintf('%sARR_TYPE', prefix);

% section
if isempty(datasection)
    datasection = lower(sprintf('.DTF_%s_sec', filter_name));
end

%% Code generation .c file
fid_c = fopen(strcat(foldername, filename_c), 'w');

%### Header, filter length and array/pointer definitions ###%
fprintf(fid_c, '#include <stdio.h>\n');
fprintf(fid_c, '#include "%s"\n\n', filename_h);
fprintf(fid_c, '#ifdef SOC_C6678\n');
fprintf(fid_c, '#include <c6x.h>\n');
fprintf(fid_c, '#endif\n');
fprintf(fid_c, '#if (%sUNIT_TEST == 1)\n', prefix);
fprintf(fid_c, '#include "%sUNIT_TEST.h"\n', prefix);
fprintf(fid_c, '#endif\n\n');

fprintf(fid_c, '/* Arrays */\n');
fprintf(fid_c, "#ifdef SOC_C6678\n");
for i = 0 : n_den-1
    fprintf(fid_c, '#pragma DATA_ALIGN(%sy%d, %d)\n',...
            prefix, i, alignment);
end
for i = 0 : n_num-1
    fprintf(fid_c, '#pragma DATA_ALIGN(%su%d, %d)\n',...
            prefix, i, alignment);
end
fprintf(fid_c, sprintf('#pragma SET_DATA_SECTION("%s")\n', datasection));
fprintf(fid_c, "#endif // SOC_C6678\n");
% DT_FILTER_ARR_TYPE DT_FILTER_y0[DT_FILTER_LEN] = {(DT_FILTER_ARR_TYPE)0.0};
tmplate_y = '%s %sy%d[%s] = {(%s)0.0};\n';
for i = 1 : n_den
    fprintf(fid_c, tmplate_y, ftype, prefix, i-1, flen, ftype);
end
tmplate_u = '%s %su%d[%s] = {(%s)0.0};\n';
for i = 1 : n_num
    fprintf(fid_c, tmplate_u, ftype, prefix, i-1, flen, ftype);
end
fprintf(fid_c, "#ifdef SOC_C6678\n");
fprintf(fid_c, '#pragma SET_DATA_SECTION()\n');
fprintf(fid_c, "#endif // SOC_C6678\n");
fprintf(fid_c, '\n\n');

fprintf(fid_c, '/* Pointers */\n');
fprintf(fid_c, "#ifdef SOC_C6678\n");
fprintf(fid_c, sprintf('#pragma SET_DATA_SECTION("%s")\n', datasection));
fprintf(fid_c, "#endif // SOC_C6678\n");
% DT_FILTER_ARR_TYPE *DT_FILTER_y0_ptr = DT_FILTER_y0
tmplate_y = '%s *%sy%d_ptr = %sy%d;\n';
for i = 1 : n_den
    fprintf(fid_c, tmplate_y, ftype, prefix, i-1, prefix, i-1);
end
tmplate_u = '%s *%su%d_ptr = %su%d;\n';
for i = 1 : n_num
    fprintf(fid_c, tmplate_u, ftype, prefix, i-1, prefix, i-1);
end
fprintf(fid_c, "#ifdef SOC_C6678\n");
fprintf(fid_c, '#pragma SET_DATA_SECTION()\n');
fprintf(fid_c, "#endif // SOC_C6678\n");
fprintf(fid_c, '\n\n');

fprintf(fid_c, '/* Coefficients */\n');
fprintf(fid_c, "#ifdef SOC_C6678\n");
fprintf(fid_c, sprintf('#pragma SET_DATA_SECTION("%s")\n', datasection));
fprintf(fid_c, "#endif // SOC_C6678\n");
% DT_FILTER_ARR_TYPE DT_FILTER_cy1 = (DT_FILTER_ARR_TYPE)0.4;
%tmplate_y = '%s %scy%d = (%s)%.16f;\n';
tmplate_y = '%s %scy%d = (%s)%.20E;\n';
for i = 2 : n_den
    fprintf(fid_c, tmplate_y, ftype, prefix, i-1, ftype, -den_f(i));
end
%tmplate_u = '%s %scu%d = (%s)%.16f;\n';
tmplate_u = '%s %scu%d = (%s)%.20E;\n';
for i = 1 : n_num
    fprintf(fid_c, tmplate_u, ftype, prefix, i-1, ftype, num_f(i));
end
fprintf(fid_c, "#ifdef SOC_C6678\n");
fprintf(fid_c, '#pragma SET_DATA_SECTION()\n');
fprintf(fid_c, "#endif // SOC_C6678\n");
fprintf(fid_c, '\n\n');


%### Functions to extract pointers ###%
% DT_FILTER_ARR_TYPE* DT_FILTER_get_y0_ptr(void)
% {
%       return DT_FILTER_y0_ptr;
% }
tmplate_y_u = '%s* %sget_%s0_ptr(void)\n{\n\treturn %s%s0_ptr;\n}\n\n\n';
fprintf(fid_c, tmplate_y_u, ftype, prefix, 'y', prefix, 'y');
fprintf(fid_c, tmplate_y_u, ftype, prefix, 'u', prefix, 'u');


%### Functions to swap pointers ###%
% void DT_FILTER_swap_y(void)
% {
% 	DT_FILTER_ARR_TYPE* tmp_y2_ptr = DT_FILTER_y2_ptr;
% 	DT_FILTER_y2_ptr = DT_FILTER_y1_ptr;
% 	DT_FILTER_y1_ptr = DT_FILTER_y0_ptr;
% 	DT_FILTER_y0_ptr = tmp_y2_ptr;
% }
tmplate_y_u_start = 'void %sswap_%s(void)\n{\n';
tmplate_y_u_swap1 = '\t%s* tmp_%s%d_ptr = %s%s%d_ptr;\n';
tmplate_y_u_swapi = '\t%s%s%d_ptr = %s%s%d_ptr;\n';
tmplate_y_u_swapn = '\t%s%s0_ptr = tmp_%s%d_ptr;\n}\n';

fprintf(fid_c, tmplate_y_u_start, prefix, 'y');
fprintf(fid_c, tmplate_y_u_swap1, ftype, 'y', n_den-1, prefix, 'y', n_den-1);
for i = n_den-1 : -1 : 1
    fprintf(fid_c, tmplate_y_u_swapi, prefix, 'y', i, prefix, 'y', i-1);
end
fprintf(fid_c, tmplate_y_u_swapn, prefix, 'y', 'y', n_den-1);
fprintf(fid_c, '\n\n');

fprintf(fid_c, tmplate_y_u_start, prefix, 'u');
fprintf(fid_c, tmplate_y_u_swap1, ftype, 'u', n_num-1, prefix, 'u', n_num-1);
for i = n_num-1 : -1 : 1
    fprintf(fid_c, tmplate_y_u_swapi, prefix, 'u', i, prefix, 'u', i-1);
end
fprintf(fid_c, tmplate_y_u_swapn, prefix, 'u', 'u', n_num-1);
fprintf(fid_c, '\n\n');

func_swap_y = sprintf('%sswap_y();', prefix);
func_swap_u = sprintf('%sswap_u();', prefix);


%### Main function ###%
% void DT_FILTER_execute(void)
% {
% 	int i;
% 
% 	DT_FILTER_swap_y();
% 
% 	for (i=0; i<DT_FILTER_LEN; i++)
% 	{
% 		DT_FILTER_y0_ptr[i] = DT_FILTER_sat(
% 				DT_FILTER_cy1 * DT_FILTER_y1_ptr[i] +
% 				//DT_FILTER_cy2 * DT_FILTER_y2_ptr[i] +//coefficient is zero
% 				DT_FILTER_cu0 * DT_FILTER_u0_ptr[i] +
% 				DT_FILTER_cu1 * DT_FILTER_u1_ptr[i] +
% 				DT_FILTER_cu2 * DT_FILTER_u2_ptr[i], DT_FILTER_MAXVAL);
% 	}
% 
% 	DT_FILTER_swap_u();
% }
func_start = sprintf('void %sexecute(void)\n{\n\tint i;\n\t\n\t%s\n\n',...
    prefix, func_swap_y);
loop_open = sprintf('\tfor (i=0; i<%s; i++)\n\t{\n\t\t%sy0_ptr[i] = %ssat(\n',...
    flen, prefix, prefix);
tmpl_loop_iter = '\t\t\t+ %sc%s%d * %s%s%d_ptr[i]\n';
tmpl_zero_loop_iter = '\t\t\t//+ %sc%s%d * %s%s%d_ptr[i]//coefficient is zero\n';
loop_close = sprintf('\t\t\t+ %scu%d * %su%d_ptr[i], %sMAXVAL);\n\t}\n\n',...
    prefix, n_num-1, prefix, n_num-1, prefix);
func_end = sprintf('\t%s\n}\n', func_swap_u);

fprintf(fid_c, func_start);
fprintf(fid_c, loop_open);
for i = 1 : n_den-1
    if den_f(i+1) == 0
        fprintf(fid_c, tmpl_zero_loop_iter, prefix, 'y', i, prefix, 'y', i);
    else
        fprintf(fid_c, tmpl_loop_iter, prefix, 'y', i, prefix, 'y', i);
    end
end
for i = 0 : n_num-2
    if num_f(i+1) == 0
        fprintf(fid_c, tmpl_zero_loop_iter, prefix, 'u', i, prefix, 'u', i);
    else
        fprintf(fid_c, tmpl_loop_iter, prefix, 'u', i, prefix, 'u', i);
    end
end
fprintf(fid_c, loop_close);
fprintf(fid_c, func_end);

% void DT_FILTER_init(void)
% {
% 	int i;
% 
% 	for (i=0; i<DT_FILTER_LEN; i++)
% 	{
% 		DT_FILTER_y0_ptr[i] = 0.0;
% 		DT_FILTER_y1_ptr[i] = 0.0;
% 		DT_FILTER_y2_ptr[i] = 0.0;
% 		DT_FILTER_u0_ptr[i] = 0.0;
% 		DT_FILTER_u1_ptr[i] = 0.0;
% 	}
% }
func_start = sprintf('\n\nvoid %sinit(void)\n{\n\tint i;\n\n', prefix);
loop_open = sprintf('\tfor (i=0; i<%s; i++)\n\t{\n',...
    flen);
loop_body_y = '\t\t%sy%d_ptr[i] = 0.0;\n';
loop_body_u = '\t\t%su%d_ptr[i] = 0.0;\n';
loop_close = sprintf('\t}\n');
func_end = sprintf('}\n');

fprintf(fid_c, func_start);
fprintf(fid_c, loop_open);
for i = 0 : n_den-1
    fprintf(fid_c, loop_body_y, prefix, i);
end
for i = 0 : n_num-1
    fprintf(fid_c, loop_body_u, prefix, i);
end
fprintf(fid_c, loop_close);
fprintf(fid_c, func_end);

% Generate unit test code - preformatted in DTF_test_code.txt
fid_t = fopen('DTF_test_code.txt');
tline = fgetl(fid_t);
fprintf(fid_c, [replace(replace(replace(tline, '%', '%%'), '\n', '\\n'), 'XXX', prefix), '\n']);
while (1)
    tline = fgetl(fid_t);
    if ~ischar(tline)
        break;
    end
    fprintf(fid_c, [replace(replace(replace(tline, '%', '%%'), '\n', '\\n'), 'XXX', prefix), '\n']);
end
fprintf(fid_t, '\n');
fclose(fid_t);
fclose(fid_c);

%% Code generation .h file
fid_h = fopen(strcat(foldername, filename_h), 'w');
ifndef_tmpl = [upper(replace(filename_h,'.','_')), '_'];
fprintf(fid_h, sprintf("#ifndef %s\n#define %s\n\n", ifndef_tmpl, ifndef_tmpl));
fprintf(fid_h, '#include "fofb_config.h"\n');
fprintf(fid_h, "#define %sUNIT_TEST    (%d)\n", prefix, print_test_data);
fprintf(fid_h, "#define %sMAXVAL    (%.10f)\n", prefix, integrator_max);
%### Top comment ###%
% /*
%  * Hard-coded vector-wise (length K=100) filter with
%  * N+1 (N=2) output and M+1 (M=2) input taps:
%  *
%  * y0 = cy1*y1+...+cyN*yN+cu0*u0+...+cuM*uM,
%  *
%  * where cyi and cui are scalar filter coefficients and
%  * yi and ui are arrays of length K.
%  *
%  */
top_com1 = sprintf('/*\n * Hard-coded vector-wise (length K=%d) filter with\n',...
    vec_len);
top_com2 = sprintf(' * N+1 (N=%d) output and M+1 (M=%d) input taps:\n',...
    n_num-1, n_den-1);
top_com3 = ' * \n * y0 = cy1*y1+...+cyN*yN+cu0*u0+...+cuM*uM,\n * \n';
top_com4 = ' * where cyi and cui are scalar filter coefficients and\n';
top_com5 = ' * yi and ui are arrays of length K.\n * \n */\n';
top_com = strcat(top_com1, top_com2, top_com3, top_com4, top_com5);
fprintf(fid_h, top_com);

fprintf(fid_h, '#define %s (%d)\n\n', flen, vec_len);
fprintf(fid_h, '#define %sXDIR (XDIR)\n\n', prefix);
fprintf(fid_h, 'typedef %s %s;\n\n\n', filter_base_type, ftype);


%### Functin get_u0_ptr with comment comment ###%
% /*
%  * DT_FILTER_ARR_TYPE* DT_FILTER_get_u0_ptr(void)
%  *
%  * Returns pointer to input vector u0. Input data needs to
%  * be written to locations 0 ... K to which get_u0_ptr is
%  * pointing to.
%  * Note: Pointer is changed after every filter call.
%  */
% DT_FILTER_ARR_TYPE* DT_FILTER_get_u0_ptr(void);
func_def = sprintf('%s* %sget_u0_ptr(void);', ftype, prefix);
get_u0_com1 = sprintf('/*\n * %s\n * \n\n', func_def);
get_u0_com2 = ' Returns pointer to input vector u0. Input data needs to\n';
get_u0_com3 = ' * be written to locations 0 ... K to which get_u0_ptr is\n';
get_u0_com4 = ' * pointing to.\n * Note: Pointer is changed after every filter call.\n */\n';
get_u0_com = strcat(get_u0_com1, get_u0_com2, get_u0_com3, get_u0_com4);
fprintf(fid_h, get_u0_com);
fprintf(fid_h, '%s\n\n\n', func_def);

func_def = sprintf('%s* %sget_y0_ptr(void);', ftype, prefix);
get_y0_com1 = sprintf('/*\n * %s\n * \n\n', func_def);
get_y0_com2 = ' Returns pointer to output vector y0. Output data needs to\n';
get_y0_com3 = ' * be read from locations 0 ... K to which get_y0_ptr is\n';
get_y0_com4 = ' * pointing to.\n * Note: Pointer is changed after every filter call.\n */\n';
get_y0_com = strcat(get_y0_com1, get_y0_com2, get_y0_com3, get_y0_com4);
fprintf(fid_h, get_y0_com);
fprintf(fid_h, '%s\n\n\n', func_def);

func_def = sprintf('void %sexecute(void);', prefix);
exec_com1 = sprintf('/*\n * void %s\n * \nExecutes filter. Should be used as follows:\n',...
    func_def);
exec_com2 = ' * 1.) Retrieve input pointer: input_data = get_u0_ptr()\n';
exec_com3 = ' * 2.) Write input data to input_data[0] ... input_data[N]\n';
exec_com4 = ' * 3.) Call DT_FILTER_execute\n';
exec_com5 = ' * 4.) Retrieve output pointer: output_data = get_y0_ptr()\n';
exec_com6 = ' */\n';
exec_com = strcat(exec_com1, exec_com2, exec_com3, exec_com4,...
    exec_com5, exec_com6);
fprintf(fid_h, exec_com);
fprintf(fid_h, '%s\n\n\n', func_def);


func_def = sprintf('void %sinit(void);', prefix);
exec_com = sprintf('/*\n * %s\n */\n', func_def);
fprintf(fid_h, exec_com);
fprintf(fid_h, func_def);
fprintf(fid_h, '\n\n');

if print_test_data == 1
    fprintf(fid_h, '\n#if (%sUNIT_TEST == 1)\n', prefix);
    fprintf(fid_h, 'int %sunit_test(void);\n', prefix);
    fprintf(fid_h, '#endif\n\n');
end

% #define DT_FILTER_min(X, Y)  ((X) < (Y) ? (X) : (Y))
% #define DT_FILTER_max(X, Y)  ((X) > (Y) ? (X) : (Y))
% #define DT_FILTER_sat(X, Y)  (DT_FILTER_min(DT_FILTER_max(X,-Y),Y))
fprintf(fid_h, sprintf('#define %smin(X, Y)  ((X) < (Y) ? (X) : (Y))\n', prefix));
fprintf(fid_h, sprintf('#define %smax(X, Y)  ((X) > (Y) ? (X) : (Y))\n', prefix));
fprintf(fid_h, sprintf('#define %ssat(X, Y)  (%smin(%smax(X,-Y),Y))\n\n', prefix, prefix, prefix));

fprintf(fid_h, sprintf("\n#endif // %s\n", ifndef_tmpl));

fclose(fid_h);


%% Print test data
if print_test_data == 1
    Ts = 1;
    nsamples = 10;
    endt = (nsamples * Ts) - Ts;
    Lsim = nsamples * Ts;
    t = 0 : Ts : endt;
    input_mat = randn(vec_len, nsamples);
    input_data = [t', input_mat'];

    % Make sure we get the expected result in simulink
    options = simset('SrcWorkspace','current');
    sim('dt_filter_initial_test_model.slx',[],options);
    
    fid = fopen(strcat(foldername, filename_test_data), 'w');
    
    ifndef_tmpl = [upper(replace(filename_test_data,'.','_')), '_'];
    fprintf(fid, sprintf("#ifndef %s\n#define %s\n\n",...
        ifndef_tmpl, ifndef_tmpl));

    fprintf(fid, 'const int %s = %d;\n',...
        sprintf('%stest_nsamples', prefix), nsamples);
    fprintf(fid, 'const int %s = %d;\n',...
        sprintf('%stest_vec_len', prefix), vec_len);

    print_dense_C_matrix(fid, input_mat', filter_base_type,...
        sprintf('%stest_input_data', prefix), true, '.gsvd_unit_test', 2);
    print_dense_C_matrix(fid, output_data(1:nsamples, :),...
        filter_base_type, sprintf('%stest_output_data', prefix), true, '.gsvd_unit_test', 2);

    fprintf(fid, sprintf("#endif // %s\n", ifndef_tmpl));
    fclose(fid);
end
end