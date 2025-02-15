function [l] = sparse_laplacian(size)
x_size = size(1);
y_size = size(2);
t_size = x_size * y_size;

sp_x_size = t_size - y_size;
sp_y_size = t_size;

sp_size = 2 * sp_x_size;
sp_row = zeros([sp_size, 1]);
sp_col = zeros([sp_size, 1]);
sp_val = zeros([sp_size, 1]);

el_count = 1;
i_row = 1;
for iy = 1:y_size
    for ix = 1:x_size - 1
        i_col_first  = ix + x_size  * (iy - 1);
        i_col_second = i_col_first + 1;

        sp_row(el_count) = i_row;
        sp_col(el_count) = i_col_first;
        sp_val(el_count) = -1;
        el_count = el_count + 1;

        sp_row(el_count) = i_row;
        sp_col(el_count) = i_col_second;
        sp_val(el_count) = 1;
        el_count = el_count + 1;

        i_row = i_row + 1;
    end
end

dx1 = sparse(sp_row, sp_col, sp_val, sp_x_size, sp_y_size);

sp_x_size = t_size - x_size;
sp_y_size = t_size;

sp_size = 2 * sp_x_size;
sp_row = zeros([sp_size, 1]);
sp_col = zeros([sp_size, 1]);
sp_val = zeros([sp_size, 1]);

el_count = 1;
i_row = 1;
for iy = 1:y_size - 1
    for ix = 1:x_size
        i_col_first = ix + (iy - 1) * x_size;
        i_col_second = i_col_first + x_size;

        sp_row(el_count) = i_row;
        sp_col(el_count) = i_col_first;
        sp_val(el_count) = -1;
        el_count = el_count + 1;

        sp_row(el_count) = i_row;
        sp_col(el_count) = i_col_second;
        sp_val(el_count) = 1;
        el_count = el_count + 1;

        i_row = i_row + 1;
    end
end

dy1 = sparse(sp_row, sp_col, sp_val, sp_x_size, sp_y_size);

l = dx1' * dx1 + dy1' * dy1;
end