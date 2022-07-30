function print_dense_C_matrix(fid, matrix, type, name, isconst, ... 
                        datasection, alignement)
    if nargin >= 6
        assert(nargin == 7);
        fprintf(fid, "#ifdef SOC_C6678\n");
        fprintf(fid, '#pragma DATA_SECTION(%s, "%s");\n',...
                name, datasection);
        fprintf(fid, '#pragma DATA_ALIGN(%s, %d);\n',...
                name, alignement);
        fprintf(fid, "#endif // SOC_C6678\n");
    end
    
    if isconst == true
        prefix = 'const ';
    else
        prefix = '';
    end
    
    [nrows, ncols] = size(matrix);
    fprintf(fid, '%s%s %s[%d] = {', prefix, type, name, (nrows * ncols));
    for row = 1 : nrows
        for col = 1 : ncols
            if strcmp(type, 'int')
                if (row == nrows) && (col == ncols)
                    fprintf(fid, '(%s)%d};\n', type, matrix(row,col));
                else
                    fprintf(fid, '(%s)%d, \n', type, matrix(row,col));
                end
            else
                if (row == nrows) && (col == ncols)
                    fprintf(fid, '(%s)%.18f};\n', type, matrix(row,col));
                else
                    fprintf(fid, '(%s)%.18f, \n', type, matrix(row,col));
                end
            end
        end
    end
end