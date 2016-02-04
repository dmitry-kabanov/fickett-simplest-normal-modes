function [znd] = compute_znd_data_at_point(x, d, q, theta, k)
    assert(isscalar(x))
    
    if x == 0.0
        u = 2*d;
        lambda = 0.0;
        u_min_d = d;
    else
        ic = 0.0;
        opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        %[~, sol] = ode45(@zndrhsfun, [0 x], ic, opts);
        [grid_, sol] = ode45(@zndrhsfun, [0 x], ic, opts);
        lambda = sol(end);
        u_min_d = sqrt(d^2 - q*lambda);
        u = d + u_min_d;
    end

    tmp = k * exp((u + q*lambda) * theta);    

    znd.u = u;
    znd.l = lambda;
    znd.w = (1 - lambda) * tmp;
    znd.dl_dx = -znd.w / d;
    znd.du_dx = -0.5 * q * znd.dl_dx / u_min_d;
    znd.dw_du = theta * znd.w;
    znd.dw_dl = tmp * ((1 - lambda) * theta * q - 1);

    
    %--------------------------------------------------------------------------
    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end