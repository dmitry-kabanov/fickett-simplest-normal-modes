function znd_sol = get_znd_sol(M, params)
q = params.q;
theta = params.theta;
d = params.d;
k = params.k;

ic = 0.0;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
znd_sol = ode45(@zndrhsfun, [0 -M], ic, opts);


    %--------------------------------------------------------------------------
    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end
