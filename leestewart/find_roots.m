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

znd_du_dx = znd_all.dl_dx(end);
znd_dl_dx = znd_all.dl_dx(end);
znd_dw_du = znd_all.dw_du(end);
znd_dw_dl = znd_all.dw_dl(end);

coeff = params.sigma / params.d;

H_real = -uprime_re * alpha_re + ...
         uprime_im * alpha_im + ...
         znd_du_dx * (alpha_re - uprime_re) + ...
         coeff * (znd_dw_du * uprime_re + ...
                  znd_dl_dx * alpha_re + ...
                  znd_dw_dl * lprime_re + ...
                  alpha_im * lprime_im - ...
                  alpha_re * lprime_re);
              
H_imag = -uprime_re * alpha_im - ...
         uprime_im * alpha_re + ...
         znd_du_dx * (alpha_im - uprime_im) + ...
         coeff * (znd_dw_du * uprime_im + ...
                  znd_dl_dx * alpha_im + ...
                  znd_dw_dl * lprime_im - ...
                  alpha_re * lprime_im - ...
                  alpha_im * lprime_re);
              
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

znd_du_dx = znd_all.dl_dx(end);
znd_dl_dx = znd_all.dl_dx(end);

coeff = params.sigma / params.d;

H_re__alpha_re = -uprime_re + znd_du_dx + coeff*(znd_dl_dx - lprime_re);
H_re__alpha_im = uprime_im + coeff*lprime_im;
H_im__alpha_re = -uprime_im - coeff*lprime_im;
H_im__alpha_im = -uprime_re + znd_du_dx + coeff*(znd_dl_dx - lprime_re);

jac = [H_re__alpha_re  H_re__alpha_im
       H_im__alpha_re  H_im__alpha_im];
end