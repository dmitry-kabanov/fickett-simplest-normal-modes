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
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
%sol = ode45(rhsfun, xspan, ic, opts);
%opts = odeset('Events', @event_check_singular_du_dx);
sol = ode45(rhsfun, xspan, ic, opts);

if sol.xe > -M
    H = 10;
else
    pert_u = sol.y(1, end);
    pert_lambda = sol.y(2, end);

    znd = compute_znd_data_at_point(-M, d, q, theta, k);

    sigma = 0.5 * q;

    r_term_1 = -znd.dw_du * pert_u;
    r_term_2 = (alpha - znd.dw_dl) * pert_lambda;
    r_term_3 = -znd.dl_dx * alpha;
    numer = r_term_1 + r_term_2 + r_term_3;
    rprime = numer / d;

    term_1 = -(alpha + znd.du_dx) * pert_u;
    term_2 = znd.du_dx * alpha;
    term_3 = -sigma * rprime;
    H = term_1 + term_2 + term_3;
end
end


%------------------------------------------------------------------------------
function rhs = rhsfun_impl(x, y, alpha, d, q, theta, k)
    u = y(1);
    lambda = y(2);
    
    rhs = zeros(2, 1);

    znd = compute_znd_data_at_point(x, d, q, theta, k);
    
    sigma = 0.5 * q;
    
    r_term_1 = -znd.dw_du * u;
    r_term_2 = (alpha - znd.dw_dl) * lambda;
    r_term_3 = -znd.dl_dx * alpha;
    numer = r_term_1 + r_term_2 + r_term_3;
    rprime = numer / d;

    u_term_1 = -(alpha + znd.du_dx) * u;
    u_term_2 = znd.du_dx * alpha;
    u_term_3 = -sigma * rprime;
    u_numer = u_term_1 + u_term_2 + u_term_3;
    rhs(1) = u_numer / (znd.u - d);
    rhs(2) = rprime;
end


%--------------------------------------------------------------------------
function [value, isterminal, direction] = event_check_singular_du_dx(~, y)
u = y(1);

value = abs(u) - 1000;
isterminal = 1;
direction = 0;
end
