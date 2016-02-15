function [x] = compute_x_at_point(lambda, params)
%COMPUTE_X_AT_POINT Compute :math:`x(\lambda)`.
f = @(s) f_impl(s, params);
x = integral(f, 0, lambda, 'RelTol', 1e-12, 'AbsTol', 1e-12);
end


%-------------------------------------------------------------------------------
function [f] = f_impl(s, params)
omega = params.k * (1 - s) * exp((u + q * lambda) * params.theta);
f = -params.d ./ omega;
end