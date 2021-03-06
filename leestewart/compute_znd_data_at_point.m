function [znd] = compute_znd_data_at_point(lambda, params)
    d = params.d;
    k = params.k;
    q = params.q;
    theta = params.theta;

    u_min_d = sqrt(d^2 - q*lambda);
    u = d + u_min_d;

    tmp = k * exp((u + q*lambda) * params.theta);

    znd.u = u;
    znd.l = lambda;
    znd.w = (1 - lambda) * tmp;

    dl_dx = -znd.w / d;
    du_dx = -0.5 * q * dl_dx / u_min_d;

    znd.du_dl = (-d / znd.w) * du_dx;

    znd.dw_du = theta * znd.w;
    znd.dw_dl = tmp * ((1 - lambda) * theta * q - 1);
end
