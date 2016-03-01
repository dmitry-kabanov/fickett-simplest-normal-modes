function [znd] = compute_znd(grid, params)
    assert(isvector(grid), 'Parameter `grid` must be a vector.')
    grid = reshape(grid, length(grid), 1);

    d = params.d;
    k = params.k;
    q = params.q;
    theta = params.theta;

    internal_grid = flipud(grid);
    ic = 0;
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [t_, sol] = ode45(@zndrhsfun, internal_grid, ic, opts);
    lambda = sol;
    assert(isvector(lambda));
    
    u_min_d = sqrt(d^2 - q*lambda);
    u = d + u_min_d;

    tmp = k * exp((u + q*lambda) * theta);    

    w =  (1 - lambda) .* tmp;
    dl_dx = -w ./ d;
    du_dx = -0.5 * q * dl_dx ./ u_min_d;
    dw_du = theta * w;
    dw_dl = tmp .* ((1 - lambda) * theta * q - 1);

    znd.u = flipud(u);
    znd.l = flipud(lambda);
    znd.du_dx = flipud(du_dx);
    znd.dl_dx = flipud(dl_dx);
    znd.w = flipud(w);
    znd.dw_du = flipud(dw_du);
    znd.dw_dl = flipud(dw_dl);
    
    
    %---------------------------------------------------------------------------
    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end
