function [znd] = compute_znd_data_at_point(x, d, q, theta, k)
    assert(isscalar(x))

    lambda = deval(znd_sol, x);
    u = d + sqrt(d^2 - q*lambda);

    tmp = k * exp((u + q*lambda) * theta);    

    w =  (1 - lambda) * tmp;
    dl_dx = -w / d;
    du_dx = 0.5 * q * w / (d * sqrt(d^2 - q*lambda));
    dw_du = theta * w;
    dw_dl = tmp * ((1 - lambda) * theta * q - 1);

    znd.u = u;
    znd.l = lambda;
    znd.du_dx = du_dx;
    znd.dl_dx = dl_dx;
    znd.dw_du = dw_du;
    znd.dw_dl = dw_dl;

    
    %--------------------------------------------------------------------------
    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end
