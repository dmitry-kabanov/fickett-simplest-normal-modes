function [znd] = compute_znd_data_at_point(lambda, params)
% Compute ZND data at point `lambda`.
d = params.d;
k = params.k;
q = params.q;
sigma = params.sigma;
theta = params.theta;

expr = d^2 - q*lambda;
assert(expr > 0)

u_min_d = sqrt(expr);
u = d + u_min_d;

tmp = k * exp((u + q*lambda) * params.theta);

znd.u = u;
znd.l = lambda;
znd.w = (1 - lambda) * tmp;

znd.dl_dx = -znd.w / d;
znd.du_dx = -sigma * znd.dl_dx / u_min_d;

znd.du_dl = (-d / znd.w) * znd.du_dx;

znd.dw_du = theta * znd.w;
znd.dw_dl = tmp * ((1 - lambda) * theta * q - 1);
end
