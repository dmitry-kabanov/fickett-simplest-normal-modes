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

[pert, sol] = compute_linearized_problem(alpha, lambda_tol, params);

if sol.xe < (1-lambda_tol)
    H = 10;
else
    pert_u = pert.u;
    pert_lambda = pert.lambda;

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



