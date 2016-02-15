function [sol] = compute_linearized_problem(alpha, grid, znd_all, params)
%COMPUTE_LINEARIZED_PROBLEM Summary of this function goes here
%   Detailed explanation goes here

rhsfun = @(t, y, idx) rhsfun_impl(t, y, idx, alpha, znd_all, params);
alpha_re = real(alpha);
alpha_im = imag(alpha);

ic = [2*alpha_re; 2*alpha_im; 0; 0];
sol = rk4(rhsfun, grid, ic);
end


%------------------------------------------------------------------------------
function rhs = rhsfun_impl(~, y, idx, alpha, znd_all, params)
    assert(length(y) == 4);
    pert_u = y(1) + 1j * y(2);
    pert_lambda = y(3) + 1j * y(4);
    pert = [pert_u; pert_lambda];
    rhs = zeros(4, 1);
    
    d = params.d;
    sigma = params.sigma;
    
    znd.u = znd_all.u(idx);
    znd.du_dx = znd_all.du_dx(idx);
    znd.dl_dx = znd_all.dl_dx(idx);
    znd.dw_du = znd_all.dw_du(idx);
    znd.dw_dl = znd_all.dw_dl(idx);

    determinant = (d^2 - znd.u * d);
    Ainv = (1 / determinant) * [-d  sigma
                                 0  znd.u - d];
                             
    I = eye(2, 2);
    C_bar = [ znd.du_dx   0
             -znd.dw_du  -znd.dw_dl];
    b = [znd.du_dx; znd.dl_dx];
    
    tmp = Ainv * (-(alpha * I + C_bar) * pert + alpha * b);

    rhs(1) = real(tmp(1));
    rhs(2) = imag(tmp(1));
    rhs(3) = real(tmp(2));
    rhs(4) = imag(tmp(2));
end


%--------------------------------------------------------------------------
function [value, isterminal, direction] = event_check_singular_du_dx(~, y)
u = y(1);

value = abs(u) - 1000;
isterminal = 1;
direction = 0;
end