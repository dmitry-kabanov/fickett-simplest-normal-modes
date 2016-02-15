function [x] = compute_x_at_point(lambda, params)
%COMPUTE_X_AT_POINT Compute :math:`x(\lambda)`.
f = @(s) f_impl(s, params);
x = integral(f, 0, lambda, 'RelTol', 1e-12, 'AbsTol', 1e-12);
end


%-------------------------------------------------------------------------------
function [f] = f_impl(s, params)
d = params.d;
q = params.q;
u = d + sqrt(d^2 - q*s);
omega = params.k * (1 - s) .* exp((u + q*s) * params.theta);
f = -d ./ omega;
end