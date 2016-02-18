function [H] = compute_boundedness_function(alpha, grid, znd_all, params)
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
sol = compute_linearized_problem(alpha, grid, znd_all, params);

pert_u_re = sol(1, end);
pert_u_im = sol(2, end);
pert_lambda_re = sol(3, end);
pert_lambda_im = sol(4, end);

pert_u = pert_u_re + 1j * pert_u_im;
pert_lambda = pert_lambda_re + 1j * pert_lambda_im;

znd_du_dx = znd_all.dl_dx(end);
znd_dl_dx = znd_all.dl_dx(end);
znd_dw_du = znd_all.dw_du(end);
znd_dw_dl = znd_all.dw_dl(end);

r_term_1 = -znd_dw_du * pert_u;
r_term_2 = (alpha - znd_dw_dl) * pert_lambda;
r_term_3 = -znd_dl_dx * alpha;
numer = r_term_1 + r_term_2 + r_term_3;
rprime = numer / params.d;

term_1 = -alpha * pert_u;
term_2 = -znd_du_dx * (pert_u - alpha);
term_3 = -params.sigma * rprime;
H = term_1 + term_2 + term_3;
end