function plot_carpet(alpha_re, alpha_im, H)
% Plot `abs(H)` in plane specified by `alpha_re` and `alpha_im`.
%
% Plot "carpet" of `abs(H)` over rectangle region of the complex plane of
% `alpha`. The range is specified by `alpha_re` and `alpha_im`.
%
% Parameters
% ----------
% alpha_re : matrix
%     Real part of `alpha`.
% alpha_im : matrix
%     Imaginary part of `alpha`.
% H : matrix
%     Complex-valued function of "boundedness condition".

[A_RE, A_IM] = meshgrid(alpha_re, alpha_im);

figure
surf(A_RE, A_IM, abs(H))
xlabel('\Re \alpha')
ylabel('\Im \alpha')
title('abs(H)')
