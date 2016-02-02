function [result] = minimize_boundedness_function(guess, M, params, znd_sol)
% Minimize boundedness function $H(\alpha)$, using `guess`.
% 
% Function $H(\alpha)$ is complex-valued, therefore, we minimize its
% absolute value.  Initial guess contains real and imaginary parts
% of $\alpha$.
% 
% Parameters
% ----------
% guess : struct with fields: [alpha_re, alpha_im]
%     Initial guess on $\alpha$.
% 
% Returns
% -------
% result : struct
%     Minimization result ??????????????????????????
alpha_re = guess.alpha_re;
alpha_im = guess.alpha_im;

% Use builtin MATLAB function for unconstrained nonlinear optimization.
fun = @(x) abs(compute_boundedness_function(x(1) + 1j * x(2), M, params, znd_sol));
opts = optimset('TolX', 1e-8, 'TolFun', 1e-8);
[x, fval, exitflag, output] = fminsearch(fun, [alpha_re, alpha_im], opts)
end
