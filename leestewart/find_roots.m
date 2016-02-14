function [ output_args ] = find_roots(guess, lambda_tol, params)
%FIND_ROOTS Find roots of boundedness condition :math:`H(\alpha)`.
%
%   To find the roots we use Newton--Raphson procedure applied to the 
%   system of two equations with two unknowns:
%
%   H_real(alpha_real, alpha_imag) = 0,
%   H_imag(alpha_real, alpha_imag) = 0
%
%   with initial guess `guess`.
%
%   Parameters
%   ----------
%   guess : array (2x1)
%       Array with two elements: initial guesses on alpha_real and
%       alpha_imag.
%
alpha = guess;
alpha_old = alpha * 1000;
tol = 1e-6;
i = 0;

while norm(alpha - alpha_old) > tol
    i = i + 1;
    alpha_old = alpha;
    f = compute_residual(alpha_old, lambda_tol, params);
    J = compute_jacobian(alpha_old, lambda_tol, params);
    alpha = alpha_old - (J \ f);
    fprintf('i = %d, |f| = %e\n', i, norm(f));
end
end


%-------------------------------------------------------------------------------
function res = compute_residual(alpha, lambda_tol, params)
% Compute residual of the system.
alpha_re = alpha(1);
alpha_im = alpha(2);

alpha_c = alpha_re + 1j*alpha_im;

[pert, ~] = compute_linearized_problem(alpha_c, lambda_tol, params);
uprime_re = real(pert.u);
uprime_im = imag(pert.u);
lprime_re = real(pert.lambda);
lprime_im = imag(pert.lambda);

znd = compute_znd_data_at_point(1-lambda_tol, params);

coeff = params.sigma / params.d;

H_real = -uprime_re * alpha_re + ...
         uprime_im * alpha_im + ...
         znd.du_dx * (alpha_re - uprime_re) + ...
         coeff * (znd.dw_du * uprime_re + ...
                  znd.dl_dx * alpha_re + ...
                  znd.dw_dl * lprime_re + ...
                  alpha_im * lprime_im - ...
                  alpha_re * lprime_re);
              
H_imag = -uprime_re * alpha_im - ...
         uprime_im * alpha_re + ...
         znd.du_dx * (alpha_im - uprime_im) + ...
         coeff * (znd.dw_du * uprime_im + ...
                  znd.dl_dx * alpha_im + ...
                  znd.dw_dl * lprime_im - ...
                  alpha_re * lprime_im - ...
                  alpha_im * lprime_re);
              
res = [H_real; H_imag];
end


%-------------------------------------------------------------------------------
function [jac] = compute_jacobian(alpha, lambda_tol, params)
% Compute jacobian of the system.
alpha_re = alpha(1);
alpha_im = alpha(2);

alpha_c = alpha_re + 1j*alpha_im;

[pert, ~] = compute_linearized_problem(alpha_c, lambda_tol, params);
uprime_re = real(pert.u);
uprime_im = imag(pert.u);
lprime_re = real(pert.lambda);
lprime_im = imag(pert.lambda);

znd = compute_znd_data_at_point(1-lambda_tol, params);

sigma = params.sigma;

H_re__alpha_re = -uprime_re + znd.du_dx + sigma*(znd.dl_dx - lprime_re);
H_re__alpha_im = uprime_im + sigma*lprime_im;
H_im__alpha_re = -uprime_im - sigma*lprime_im;
H_im__alpha_im = -uprime_re + znd.du_dx + sigma*(znd.dl_dx - lprime_re);

jac = [H_re__alpha_re  H_re__alpha_im
       H_im__alpha_re  H_im__alpha_im];
end