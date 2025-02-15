function [k] = wendlandC2(r)
% The Wendland Kernel that lives in a second order Sobelov space
s = 1; % support
k = zeros(size(r)); % values outside support
k(r < s) = (1 + 2 * r(r < s).^2) .* sqrt(1 - r(r < s).^2) + ...
    3 * r(r < s).^2 .* (log(r(r < s)) - ...
    log(1 + sqrt(1 - r(r < s).^2))); % values within support
end