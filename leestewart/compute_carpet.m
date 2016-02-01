function [ALPHA_RE, ALPHA_IM, H] = compute_carpet(params, M, cp)
%COMPUTE_CARPET     Compute a "carpet" of "boundedness condition".
%
%   To do linear stability analysis of the Fickett--Faria model by the
%   Lee--Stewart approach [1], [2], it takes computation of a "carpet".
%   That is, computation of the function `H` of the "boundedness condition"
%   on the grid of complex variable `alpha`. Looking at the plot of
%   `abs(H)`, it is possible to find local minima, which give initial
%   guesses for `alpha`, which can be used to find roots of `H` to the
%   machine precision.
%
%   Parameters
%   ----------
%   params : struct
%       Free and dependent parameters of the ZND solution.
%       q : float
%           Heat release.
%       theta : float
%           Activation energy.
%       d : float
%           Detonation velocity
%       k : float
%           Arrhenius prefactor.
%   M : float
%       The length of the reaction zone. Domain will be [-M 0].
%   cp : struct
%       Carpet parameters.
%       n_re : int
%           Resolution of the carpet along real axis.
%       n_im : int
%           Resolution of the carpet along imaginary axis.
%       lb_re, ub_re : float
%           Lower and upper bounds on the real axis.
%       lb_im, ub_im : float
%           Lower and upper bounds on the imaginary axis.
%
%   Returns
%   -------
%   ALPHA_RE : array (`N_RE` by `N_IM`)
%       Real part of `alpha` grid.
%   ALPHA_IM : array (`N_RE` by `N_IM`)
%       Imaginary part of `alpha` grid.
%   H : array (`N_RE` by `N_IM`)
%       Values of the "boundedness-condition" function `H` on the grid.
%   
%   References
%   ----------
%   [1] Lee, H, I. and Stewart, D. S.
%       Calculation of linear detonation instability: one-dimensional
%       instability of plane detonation.
%       Journal of Fluid Mechanics, vol. 216, pp. 103--132, 1990.
%       doi: 10.1017/S0022112090000362
%   [2] Faria, Luiz M. and Kasimov, Aslan R. and Rosales, Rodolfo R.
%       Theory of Weakly Nonlinear Self-sustained Detonations.
%       Journal of Fluid Mechanics, vol. 784, pp. 163--198, 2015.
%       doi: 10.1017/jfm.2015.577

alpha_re = linspace(cp.lb_re, cp.ub_re, cp.n_re);
alpha_im = linspace(cp.lb_im, cp.ub_im, cp.n_im);

H = zeros(cp.n_re, cp.n_im);

for j = 1:cp.n_im
    parfor i = 1:cp.n_re
        alpha = alpha_re(i) + 1j * alpha_im(j);
        
        H(i, j) = compute_boundedness_function(alpha, M, params);
        fprintf('%d, %d\n', i, j);
    end
end

[ALPHA_RE, ALPHA_IM] = meshgrid(alpha_re, alpha_im);
end
