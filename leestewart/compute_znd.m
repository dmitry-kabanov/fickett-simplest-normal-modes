function [znd] = compute_znd(grid, params)
    d = params.d;
    k = params.k;
    q = params.q;
    theta = params.theta;
    
    ic = 0.0;
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [t_, sol] = ode45(@zndrhsfun, grid, ic, opts);
    lambda = sol;
    assert(isvector(lambda));
    
    u_min_d = sqrt(d^2 - q*lambda);
    u = d + u_min_d;

    tmp = k * exp((u + q*lambda) * theta);    

    w =  (1 - lambda) .* tmp;
    dl_dx = -w ./ d;
    du_dx = -0.5 * q * dl_dx ./ u_min_d;

    znd.u = u;
    znd.l = lambda;
    znd.du_dx = du_dx;
    znd.dl_dx = dl_dx;
    znd.w = w;
    znd.dw_du = theta * znd.w;
    znd.dw_dl = tmp .* ((1 - lambda) * theta * q - 1);
    
    
    %---------------------------------------------------------------------------
    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end
