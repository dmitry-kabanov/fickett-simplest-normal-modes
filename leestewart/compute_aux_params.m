function [params] = compute_aux_params(q, theta)
% Compute auxiliary parameters of ZND solution using free parameters.
% 
% Auxiliary parameters are ZND detonation velocity and Arrhenius prefactor,
% which can be defined from free parameters `q` and `theta`.
% 
% Parameters
% ----------
% q : float
%     Heat release.
% theta : float
%     Activation energy.
% 
% Returns
% -------
% params : struct
%     Structure with free and auxiliary parameters defining ZND solution.

d = sqrt(q);

kfun = @(x) kfun_impl(x, d, q, theta);
k = integral(kfun, 0, 0.5);
params.q = q;
params.theta = theta;
params.d = d;
params.k = k;
end


%------------------------------------------------------------------------------
function res = kfun_impl(x_, d, q, theta)
    u_tmp = d + sqrt(d^2 - q*x_);
    denom = (1 - x_) .* exp((u_tmp + q*x_) * theta);
    res = d ./ denom;
end
