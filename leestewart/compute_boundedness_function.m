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
pert_l_re = sol(3, end);
pert_l_im = sol(4, end);

pert_u = pert_u_re + 1j * pert_u_im;
pert_l = pert_l_re + 1j * pert_l_im;

znd_u = znd_all.u(end);
znd_dw_dl = znd_all.dw_dl(end);

sigma = params.sigma;

term_1 = alpha * (znd_u * pert_u + sigma * pert_l);
term_2 = -sigma * znd_dw_dl * pert_l;
H = term_1 + term_2;
end