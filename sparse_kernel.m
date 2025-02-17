function [G] = sparse_kernel(obs_X, obs_Y, grid_X_vector, grid_Y_vector, fill_space)
% Get bounding box on geodetic coorinates
obs_X_min = min(obs_X, [], "all");
obs_X_max = max(obs_X, [], "all");
obs_Y_min = min(obs_Y, [], "all");
obs_Y_max = max(obs_Y, [], "all");

% Get bounding box on geodetic coordinates
grid_X_min = min(grid_X_vector, [], "all");
grid_X_max = max(grid_X_vector, [], "all");
grid_Y_min = min(grid_Y_vector, [], "all");
grid_Y_max = max(grid_Y_vector, [], "all");

% Geodetic X & Y => [0, 1]x[0, 1]
X_min = min(grid_X_min, obs_X_min);
X_max = max(grid_X_max, obs_X_max);
Y_min = min(grid_Y_min, obs_Y_min);
Y_max = max(grid_Y_max, obs_Y_max);

X_dist = X_max - X_min;
Y_dist = Y_max - Y_min;
rsc_dist = max(X_dist, Y_dist);

rsc_obs_X = (obs_X - X_min) / rsc_dist;
rsc_obs_Y = (obs_Y - Y_min) / rsc_dist;

rsc_grid_X_vector = (grid_X_vector - X_min) / rsc_dist;
rsc_grid_Y_vector = (grid_Y_vector - Y_min) / rsc_dist;

rsc_grid_X_step = rsc_grid_X_vector(2) - rsc_grid_X_vector(1);
rsc_grid_Y_step = rsc_grid_Y_vector(2) - rsc_grid_Y_vector(1);

rsc_fill_space = fill_space / rsc_dist;

rsc_X_min = min(rsc_grid_X_vector, [], "all");
rsc_Y_min = min(rsc_grid_Y_vector, [], "all");

sp_row = [];
sp_col = [];
sp_val = [];

obs_size = numel(obs_X);
grid_X_dim_size = numel(grid_X_vector);
grid_Y_dim_size = numel(grid_Y_vector);
grid_size = grid_X_dim_size * grid_Y_dim_size;
ix_step = floor(rsc_fill_space / rsc_grid_X_step);
iy_step = floor(rsc_fill_space / rsc_grid_Y_step);

if (ix_step == 0)
    warning("The radius of the kernel is less than a x grid cell spacing!")
end
if (iy_step == 0)
    warning("The radius of ther kenel is less than a y grid cell spacing!")
end

for k = 1:obs_size
	rsc_obs_X_k = rsc_obs_X(k);
	rsc_obs_Y_k = rsc_obs_Y(k);

	rsc_obs_dx = rsc_obs_X_k - rsc_X_min;
	ix_middle = floor(rsc_obs_dx / rsc_grid_X_step) + 1;

	ix_lower = ix_middle - ix_step;
	ix_upper = ix_middle + ix_step;

	rsc_obs_dy = rsc_obs_Y_k - rsc_Y_min;
	iy_middle = floor(rsc_obs_dy / rsc_grid_Y_step) + 1;

	iy_lower = iy_middle - iy_step;
	iy_upper = iy_middle + iy_step;

	if iy_lower < 1
        if iy_upper < 1
            continue
        end
		iy_lower = 1;
	end
	if iy_upper > grid_Y_dim_size
        if iy_lower > grid_Y_dim_size
            continue
        end
		iy_upper = grid_Y_dim_size;
    end
    if iy_lower > iy_upper
        continue
    end
	if ix_lower < 1
        if ix_upper < 1
            continue
        end
		ix_lower = 1;
	end
	if ix_upper > grid_X_dim_size
        if ix_lower > grid_X_dim_size
            continue
        end
	    ix_upper = grid_X_dim_size;
    end
    if ix_lower > ix_upper
        continue
    end

	for ix_k = ix_lower:ix_upper
		rsc_grid_X_vector_k = rsc_grid_X_vector(ix_k);
		rsc_grid_dx = rsc_grid_X_vector_k - rsc_obs_X_k;
		
        for iy_k = iy_lower:iy_upper
			rsc_grid_Y_vector_k = rsc_grid_Y_vector(iy_k);
			rsc_grid_dy = rsc_grid_Y_vector_k - rsc_obs_Y_k;

			rsc_grid_d = sqrt(rsc_grid_dx^2 + rsc_grid_dy^2);
			ker = wendlandC2(rsc_grid_d / rsc_fill_space);

			if ker > 0
				sp_row(end + 1) = k;
				sp_col(end + 1) = ix_k + grid_X_dim_size * (iy_k - 1);
				sp_val(end + 1) = ker;
			end
        end
    end
end
G = sparse(sp_row, sp_col, sp_val, obs_size, grid_size);
end