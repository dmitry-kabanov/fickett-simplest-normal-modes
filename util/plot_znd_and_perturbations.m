function plot_znd_and_perturbations(grid, znd_all, pert)
%PLOT_ZND_AND_PERTURBATIONS Plot ZND solution and perturbations.
%
% ZND solution, real parts and imaginary parts of perturbations are plotted
% vertically stacked, so that it is convenient to see that perturbations
% vanish as ZND solution approached equilibrium.
%
% Parameters
% ----------
% grid: array (Nx1)
%     Grid
% znd_all: struct
%     Structure containing ZND solution on the `grid`.
% pert: array (4xN)
%     Array containing real and imaginary parts of perturbations in this
%     order: real(u), imag(u), real(lambda), imag(lambda).
%
figure
subplot(3, 1, 1, 'align');
hold on
plot(grid, znd_all.u, '-')
plot(grid, znd_all.l, '--')
plot(grid, znd_all.w, '-.')
ylabel('ZND solution')
legend('u', '\lambda', '\omega', 'Location', 'northwest')

subplot(3, 1, 2, 'align');
hold on
plot(grid, pert(1, :))
plot(grid, pert(3, :))
ylabel('Real part')
legend('Re u', 'Re \lambda', 'Location', 'southwest');

subplot(3, 1, 3, 'align');
hold on
plot(grid, pert(2, :))
plot(grid, pert(4, :))
xlabel('x')
ylabel('Imaginary part')
legend('Im u', 'Im \lambda', 'Location', 'southwest');
end