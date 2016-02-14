function [x] = compute_x_at_point(lambda, params)
%COMPUTE_X_AT_POINT Compute :math:`x(\lambda)`.
f = @(s) f_impl(s, params);
x = integral(f, 0, lambda, 'RelTol', 1e-12, 'AbsTol', 1e-12);
end


%-------------------------------------------------------------------------------
function [f] = f_impl(s, params)
znd = compute_znd_data_at_point(s, params);
omega = znd.w;
f = -params.d ./ omega;
end

