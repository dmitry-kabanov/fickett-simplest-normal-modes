% Check that we compute ZND structure correctly when we use function
% `compute_znd_data`.
close all;
clear;

% Physical free parameters.
q = 1.7; theta = 2.4;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

d = params.d;
k = params.k;

% Domain length.
M = 13;

x = linspace(-M, 0, 1000);
u = zeros('like', x);
lambda = zeros('like', x);
du_dx = zeros('like', x);
dl_dx = zeros('like', x);

for i = 1:length(x)
    znd = compute_znd_data_at_point(x(i), d, q, theta, k);
    u(i) = znd.u;
    lambda(i) = znd.l;
    du_dx(i) = znd.du_dx;
    dl_dx(i) = znd.dl_dx;
end
