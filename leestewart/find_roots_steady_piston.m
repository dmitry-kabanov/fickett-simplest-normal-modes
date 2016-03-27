function [result] = find_roots_steady_piston(guess, grid, znd_all, params, tol)
% Find eigenvalue using shooting method with steady-piston boundary condition.
% 
% Steady-piston boundary condition is :math:`u' = 0`, where :math:`u'` is a
% complex-valued function, therefore, we minimize its absolute value.
% Initial guess contains real and imaginary parts of $\alpha$.
% 
% Parameters
% ----------
% guess : struct with fields: [alpha_re, alpha_im]
%     Initial guess on $\alpha$.
% grid : array
%     Grid.
% znd_all : struct
%     Structure containing ZND solution and auxiliary quantities computed
%     on the `grid`.
% params : struct
%     Parameters of the problem.
% 
% Returns
% -------
% res : struct
%     Minimization result. Fields:
%         root : array (1x2)
%             Array containing minimizing values of :math:`\alpha`.
%         fval : float
%             Minimum value of :math:`\abs(u')`.
%         output : struct
%             Structure containing information about the optimization.
%             Provided by `fminsearch`.
% Use builtin MATLAB function for unconstrained nonlinear optimization.
fun = @(x) fun_impl(x, grid, znd_all, params);
opts = optimoptions('fsolve', 'TolX', tol, 'TolFun', tol, 'Display', 'iter');
[x, fval, exitflag, output] = fsolve(fun, guess, opts);

%assert(exitflag == 1, '`fminsearch` terminated with not successful result');

result.root = x;
result.fval = fval;
result.exitflag = exitflag;
result.output = output;
end


%-------------------------------------------------------------------------------
function res = fun_impl(given, grid, znd_all, params)
% Compute residual of the system.
ic = [given(1); given(2); given(3); given(4)];

alpha_re = given(5);
alpha_im = given(6);

alpha_c = alpha_re + 1j*alpha_im;

sol = compute_linearized_problem(alpha_c, ic, grid, znd_all, params);
uprime_re = sol(1, end);
uprime_im = sol(2, end);
lprime_re = sol(3, end);
lprime_im = sol(4, end);

rhcond = [uprime_re - 2*alpha_re;
          uprime_im - 2*alpha_im;
          lprime_re;
          lprime_im];

res = norm(rhcond, 2);
end