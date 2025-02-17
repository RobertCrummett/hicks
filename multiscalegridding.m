%% Hicks Dome Gridding
% R Nate Crummett, 2/15/25
 clear

%% Load observed (obs) data
obs_path = "combined_data_rad_mag.csv";
obs_xyz = readtable(obs_path, "ReadVariableNames", false);
obs_xyz.Properties.VariableNames = {'X', 'Y', 'K', 'T', 'U', 'TMI'};

obs_X = obs_xyz.X;
obs_Y = obs_xyz.Y;
obs_tmi = obs_xyz.TMI;

% Remove tie lines
ties = 4200;
obs_X = obs_X(1:end-ties);
obs_Y = obs_Y(1:end-ties);
obs_tmi = obs_tmi(1:end-ties);

% Get bounding box on geodetic coorinates
obs_X_min = min(obs_X, [], "all");
obs_X_max = max(obs_X, [], "all");
obs_Y_min = min(obs_Y, [], "all");
obs_Y_max = max(obs_Y, [], "all");

% Remove DC component
obs_tmi_mean = mean(obs_tmi);
obs_tfa = obs_tmi - obs_tmi_mean;

obs_size = length(obs_X);

%% Load model (grid) data
grid_path = "mag_grid_data.csv";
grid_xyz = readtable(grid_path, "ReadVariableNames", false);
grid_xyz.Properties.VariableNames = {'X', 'Y', 'TMI'};

grid_X_dim_size = 500;
grid_Y_dim_size = 500;

grid_shape = [grid_X_dim_size, grid_Y_dim_size];
grid_size = prod(grid_shape);

grid_X = reshape(grid_xyz.X, grid_shape);
grid_Y = reshape(grid_xyz.Y, grid_shape);

grid_X_vector = grid_X(:,1);
grid_Y_vector = grid_Y(1,:)';

% % Get bounding box on geodetic coordinates
% grid_X_min = min(grid_X_vector, [], "all");
% grid_X_max = max(grid_X_vector, [], "all");
% grid_Y_min = min(grid_Y_vector, [], "all");
% grid_Y_max = max(grid_Y_vector, [], "all");

% grid_X_step = grid_X_vector(2) - grid_X_vector(1);
% grid_Y_step = grid_Y_vector(2) - grid_Y_vector(1);

%% Rescale (rsc) observations and grid locations
% Geodetic X & Y => [0, 1]x[0, 1]
% X_min = min(grid_X_min, obs_X_min);
% X_max = max(grid_X_max, obs_X_max);
% Y_min = min(grid_Y_min, obs_Y_min);
% Y_max = max(grid_Y_max, obs_Y_max);
% 
% X_dist = X_max - X_min;
% Y_dist = Y_max - Y_min;
% rsc_dist = max(X_dist, Y_dist);
% 
% rsc_obs_X = (obs_X - X_min) / rsc_dist;
% rsc_obs_Y = (obs_Y - Y_min) / rsc_dist;
% 
% rsc_grid_X_vector = (grid_X_vector - X_min) / rsc_dist;
% rsc_grid_Y_vector = (grid_Y_vector - Y_min) / rsc_dist;
% 
% rsc_grid_X_step = rsc_grid_X_vector(2) - rsc_grid_X_vector(1);
% rsc_grid_Y_step = rsc_grid_Y_vector(2) - rsc_grid_Y_vector(1);
% 
% % Get bounding box in [0,1]x[0,1]
% rsc_X_min = min(rsc_grid_X_vector, [], "all");
% rsc_X_max = max(rsc_grid_X_vector, [], "all");
% rsc_Y_min = min(rsc_grid_Y_vector, [], "all");
% rsc_Y_max = max(rsc_grid_Y_vector, [], "all");

%% Scale 1 - Dome scale ~ 20 km
scale1_percent = 0.1;

scale1_choices = randsample(1:obs_size, int32(scale1_percent * obs_size), false);
scale1_mask = zeros(obs_size, 1, "logical");
for s1c = scale1_choices
    scale1_mask(s1c) = true;
end

scale1_obs_X = obs_X(scale1_mask);
scale1_obs_Y = obs_Y(scale1_mask);
scale1_tfa = obs_tfa(scale1_mask);

%figure()
%scatter(scale1_X, scale1_Y, 10, scale1_tfa, "filled")

%% Scale 1 - Division
scale1_div_X = 20;
scale1_div_Y = 20;

scale1_X_dim_size = floor(grid_X_dim_size / scale1_div_X);
scale1_Y_dim_size = floor(grid_Y_dim_size / scale1_div_Y);

scale1_grid_tfa = zeros(grid_shape);

if scale1_X_dim_size == 0
    error("The sub grid x dimension is 0")
end
if scale1_Y_dim_size == 0
    error("The sub grid y dimension is 0")
end

for i=1:scale1_div_X
    scale1_div_X_start = (i-1) * scale1_X_dim_size + 1;
    scale1_div_X_end = i * scale1_X_dim_size;
    if i == scale1_div_X
        scale1_div_X_end = grid_X_dim_size;
    end

    for j=1:scale1_div_Y
        scale1_div_Y_start = (j-1) * scale1_Y_dim_size + 1;
        scale1_div_Y_end = j * scale1_Y_dim_size;
        if j == scale1_div_Y
            scale1_div_Y_end = grid_Y_dim_size;
        end

        scale1_div_shape = [scale1_div_X_end-scale1_div_X_start+1,...
            scale1_div_Y_end-scale1_div_Y_start+1];
        scale1_div_size = prod(scale1_div_shape);

        % figure()
        % scatter(scale1_obs_X, scale1_obs_Y, 10, scale1_tfa)
        % hold on
        % scatter(grid_X_vector(scale1_div_X_start:scale1_div_X_end), ...
        %     grid_Y_vector(scale1_div_Y_start:scale1_div_Y_end), 10, "b")
        % hold off
        % uiwait()

        scale1_fill = 4e3;
        scale1_G = sparse_kernel(scale1_obs_X, scale1_obs_Y, ...
            grid_X_vector(scale1_div_X_start:scale1_div_X_end), ...
            grid_Y_vector(scale1_div_Y_start:scale1_div_Y_end), scale1_fill);
        
        scale1_L = sparse_laplacian(scale1_div_shape);
        scale1_alpha = 1e3;
        scale1_mof = speye(scale1_div_size) + scale1_alpha * scale1_L;

        scale1_beta = 1e1;
        scale1_div_tfa = (scale1_G' * scale1_G + scale1_beta * scale1_mof) \ (scale1_G' * scale1_tfa);
        scale1_div_tfa = reshape(scale1_div_tfa, scale1_div_shape);

        scale1_grid_tfa( ...
            scale1_div_Y_start:scale1_div_Y_end, ...
            scale1_div_X_start:scale1_div_X_end) = scale1_div_tfa;

        fprintf("%d, %d\n", i, j);
    end
end

if sum(scale1_tfa, "all") == grid_size
    disp("Yes!")
end
%% Scale 1 - Inversion

% scale1_fill = 2e3;
% scale1_G = sparse_kernel(scale1_obs_X, scale1_obs_Y, grid_X_vector, grid_Y_vector, scale1_fill);
% 
% scale1_L = sparse_laplacian(grid_shape);
% scale1_alpha = 1e3;
% scale1_mof = speye(grid_size) + scale1_alpha * scale1_L;
% 
% scale1_bata = 1e-3;
% scale1_tfa = (scale1_G' * scale1_G + scale1_beta * scale1_mof) \ (scale1_G' * scale1_tfa);
% 
% scale1_tfa = reshape(scale1_tfa, grid_shape);