function [res] = minimize_boundedness_function(guess, grid, znd_all, params)
% Minimize boundedness function :math:`H(\alpha)`, using `guess`.
% 
% Function $H(\alpha)$ is complex-valued, therefore, we minimize its
% absolute value.  Initial guess contains real and imaginary parts
% of $\alpha$.
% 
% Parameters
% ----------
% guess : struct with fields: [alpha_re, alpha_im]
%     Initial guess on $\alpha$.
% lambda_tol : float
%     Domain length.
% params : struct
%     Parameters of the problem.
% 
% Returns
% -------
% res : struct
%     Minimization result ??????????????????????????
alpha_re = guess.alpha_re;
alpha_im = guess.alpha_im;

% Use builtin MATLAB function for unconstrained nonlinear optimization.
fun = @(x) compute_objective(x, grid, znd_all, params);
opts = optimoptions('fsolve', ...
                    'TolX', 1e-15, ...
                    'TolFun', 1e-14, ...
                    'Algorithm', 'trust-region-reflective', ...
                    'FinDiffType', 'central', ...
                    'Display', 'iter');

problem.objective = fun;
problem.x0 = [alpha_re, alpha_im];
problem.solver = 'fsolve';
problem.options = opts;
[x, fval, exitflag, output] = fsolve(problem);

res.x = x;
res.fval = fval;
res.exitflag = exitflag;
res.output = output;
end

%-------------------------------------------------------------------------------
function residual = compute_objective(x, grid, znd_all, params)
alpha_c = x(1) + 1j * x(2);
H = compute_boundedness_function(alpha_c, grid, znd_all, params);

residual = [real(H); imag(H)];
end