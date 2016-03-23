function [root, norm_f] = find_roots(guess, grid, znd_all, params, tol)
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
i = 0;

while norm(alpha - alpha_old) > tol
    i = i + 1;
    alpha_old = alpha;
    f = compute_residual(alpha_old, grid, znd_all, params);
    J = compute_jacobian(alpha_old, grid, znd_all, params);
    alpha = alpha_old - (J \ f);
    fprintf('i = %d, |f| = %e\n', i, norm(f));
end

root = alpha;
norm_f = norm(f);
end


%-------------------------------------------------------------------------------
function res = compute_residual(alpha_arr, grid, znd_all, params)
% Compute residual of the system.
alpha_re = alpha_arr(1);
alpha_im = alpha_arr(2);

alpha_c = alpha_re + 1j*alpha_im;

sol = compute_linearized_problem(alpha_c, grid, znd_all, params);
uprime_re = sol(1, end);
uprime_im = sol(2, end);
lprime_re = sol(3, end);
lprime_im = sol(4, end);

znd_u = znd_all.u(end);
znd_dw_dl = znd_all.dw_dl(end);

sigma = params.sigma;

H_real = -znd_u * uprime_im * alpha_im + ...
          znd_u * uprime_re * alpha_re - ...
          sigma * alpha_im * lprime_im - ...
          sigma * znd_dw_dl * lprime_re + ...
          sigma * alpha_re * lprime_re;
              
H_imag = znd_u * uprime_re * alpha_im + ...
         znd_u * uprime_im * alpha_re - ...
         sigma * znd_dw_dl * lprime_im + ...
         sigma * alpha_re * lprime_im + ...
         sigma * alpha_im * lprime_re;
              
res = [H_real; H_imag];
end


%-------------------------------------------------------------------------------
function [jac] = compute_jacobian(alpha_arr, grid, znd_all, params)
% Compute jacobian of the system.
alpha_re = alpha_arr(1);
alpha_im = alpha_arr(2);

alpha_c = alpha_re + 1j*alpha_im;

sol = compute_linearized_problem(alpha_c, grid, znd_all, params);
uprime_re = sol(1, end);
uprime_im = sol(2, end);
lprime_re = sol(3, end);
lprime_im = sol(4, end);

znd_u = znd_all.u(end);

sigma = params.sigma;

H_re__alpha_re = znd_u * uprime_re + sigma * lprime_re;
H_re__alpha_im = -znd_u * uprime_im - sigma * lprime_im;
H_im__alpha_re = znd_u * uprime_im + sigma * lprime_im;
H_im__alpha_im = znd_u * uprime_re + sigma * lprime_re;

jac = [H_re__alpha_re  H_re__alpha_im
       H_im__alpha_re  H_im__alpha_im];
end