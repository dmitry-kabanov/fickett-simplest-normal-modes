function [H] = compute_boundedness_function(alpha, lambda_tol, params)
% Compute "boundedness function" :math:`H(\alpha)`.
% 
% Parameters
% ----------
% alpha : complex
%     Growth rate and frequency of perturbation.
% lambda_tol : float
%     Tolerance of :math:`lambda`, that is, closeness of the reaction
%     progress variable to unity in the end of the reaction zone.
% params : struct
%     Parameters of the problem and the ZND solution.
% 
% Returns
% -------
% H : complex
%     Value of the "boundedness function".

ic = [2*alpha 0];
xspan = [0 1-lambda_tol];

rhsfun = @(x, y) rhsfun_impl(x, y, alpha, params);
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, ...
              'Events', @event_check_singular_du_dx);
sol = ode15s(rhsfun, xspan, ic, opts);

if sol.xe < (1-lambda_tol)
    H = 10;
else
    pert_u = sol.y(1, end);
    pert_lambda = sol.y(2, end);

    znd = compute_znd_data_at_point(1-lambda_tol, params);

    r_term_1 = -znd.dw_du * pert_u;
    r_term_2 = (alpha - znd.dw_dl) * pert_lambda;
    r_term_3 = -alpha;
    numer = r_term_1 + r_term_2 + r_term_3;
    rprime = numer / params.d;

    term_1 = -(alpha + znd.du_dl) * pert_u;
    term_2 = znd.du_dl * alpha;
    term_3 = params.d * params.sigma * rprime / znd.w;
    H = term_1 + term_2 + term_3;
end
end


%------------------------------------------------------------------------------
function rhs = rhsfun_impl(x, y, alpha, params)
    u = y(1);
    lambda = y(2);
    rhs = zeros(2, 1);

    d = params.d;
    znd = compute_znd_data_at_point(x, params);
    
    r_term_1 = -znd.dw_du * u;
    r_term_2 = (alpha - znd.dw_dl) * lambda;
    r_term_3 = -alpha;
    numer = r_term_1 + r_term_2 + r_term_3;
    rprime = numer / d;

    u_term_1 = -(alpha + znd.du_dl) * u;
    u_term_2 = znd.du_dl * alpha;
    u_term_3 = d * params.sigma * rprime / znd.w;
    u_numer = u_term_1 + u_term_2 + u_term_3;
    rhs(1) = u_numer / (znd.u - d);
    rhs(2) = rprime;

    rhs = (-d / znd.w) * rhs;
end


%--------------------------------------------------------------------------
function [value, isterminal, direction] = event_check_singular_du_dx(~, y)
u = y(1);

value = abs(u) - 1000;
isterminal = 1;
direction = 0;
end
