function [H] = compute_boundedness_function(alpha, M, params)
% Compute "boundedness function" $H(\alpha)$.
% 
% Parameters
% ----------
% alpha : complex
%     Growth rate and frequency of perturbation.
% M : float
%     Domain length.
% params : struct
%     Parameters of the ZND solution.
% 
% Returns
% -------
% H : complex
%     Value of the "boundedness function".
q = params.q;
theta = params.theta;
d = params.d;
k = params.k;

ic = [2*alpha 0];
xspan = [0 -M];

rhsfun = @(x, y) rhsfun_impl(x, y, alpha, d, q, theta, k);
[~, y] = ode45(rhsfun, xspan, ic);
pert_u = y(1);
pert_lambda = y(2);

znd = compute_znd_data(-M, d, q, theta, k);

sigma = 0.5 * q;
term_1 = -(alpha + znd.du_dx + sigma * znd.dw_du) * pert_u;
term_2 = sigma * znd.dw_dl * pert_lambda;
term_3 = alpha * znd.du_dx;
H = term_1 + term_2 + term_3;
end


%------------------------------------------------------------------------------
function rhs = rhsfun_impl(x, y, alpha, d, q, theta, k)
    u = y(1);
    lambda = y(2);
    
    rhs = zeros(2, 1);

    znd = compute_znd_data(x, d, q, theta, k);
    
    sigma = 0.5 * q;
    
    term_1 = -(alpha + znd.du_dx + sigma * znd.dw_du) * u;
    term_2 = sigma * znd.dw_dl * lambda;
    term_3 = alpha * znd.du_dx;
    numer = term_1 + term_2 + term_3;
    denom = znd.u - d;
    rhs(1) = numer / denom;

    term_1 = -znd.dw_du * u;
    term_2 = (alpha - znd.dw_dl) * lambda;
    term_3 = - alpha * znd.dl_dx;
    rhs(2) = term_1 + term_2 + term_3;
end


%------------------------------------------------------------------------------
function [znd] = compute_znd_data(x, d, q, theta, k)
    if x == 0.0
        u = 2*d;
        lambda = 0.0;
    else
        ic = 0.0;
        [~, sol] = ode45(@zndrhsfun, [0 x], ic);
        lambda = sol(end);
        u = d + sqrt(d^2 - q*lambda);
    end

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


    function res = zndrhsfun(~, y)
       lambda_ = y(1);
       u_ = d + sqrt(d^2 - q*lambda_);
       omega_ = k * (1 - lambda_) * exp((u_ + q*lambda_) * theta);
       res = -omega_ / d;
    end
end
