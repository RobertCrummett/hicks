function [dx, dy] = mag_derivatives(tfa)
% Lanczos least squares derivatives in x and y direction
shape = size(tfa);
dx = zeros(shape);
dy = zeros(shape);

s0 = tfa(1:end-4,:);
s1 = tfa(2:end-3,:);
s3 = tfa(4:end-1,:);
s4 = tfa(5:end,:);

dx(3:end-2,:) = (-2 * s0 - s1 + s3 + 2 * s4) / 10;
dx(1,:) = dx(3,:);
dx(2,:) = dx(3,:);
dx(end-1,:) = dx(end-2,:);
dx(end,:) = dx(end-2,:);

s0 = tfa(:,1:end-4);
s1 = tfa(:,2:end-3);
s3 = tfa(:,4:end-1);
s4 = tfa(:,5:end);

dy(:,3:end-2) = (-2 * s0 - s1 + s3 + 2 * s4) / 10;
dy(:,1) = dy(:,3);
dy(:,2) = dy(:,3);
dy(:,end-1) = dy(:,end-2);
dy(:,end) = dy(:,end-2);
end