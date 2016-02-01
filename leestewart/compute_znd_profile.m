function [grid, znd] = compute_znd_profile(M, params)
% Compute ZND profile.

q = params.q;
theta = params.theta;
d = params.d;
k = params.k;

ic = 0.0;
opts = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
[grid, sol] = ode45(@zndrhsfun, [0 -M], ic, opts);
lambda = sol;
u = d + sqrt(d^2 - q*lambda);

tmp = k * exp((u + q*lambda) * theta);    

w =  (1 - lambda) .* tmp;
dl_dx = -w / d;
du_dx = 0.5 * q * w ./ (d * sqrt(d^2 - q*lambda));

znd.u = u;
znd.l = lambda;
znd.du_dx = du_dx;
znd.dl_dx = dl_dx;


    %--------------------------------------------------------------------------
    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end
