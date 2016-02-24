function plot_znd_and_perturbations(matfile, picfile)
%PLOT_ZND_AND_PERTURBATIONS Plot ZND solution and perturbations.
%
% ZND solution, real parts and imaginary parts of perturbations are plotted
% vertically stacked, so that it is convenient to see that perturbations
% vanish as ZND solution approached equilibrium.
%
% Parameters
% ----------
% matfile : str
%     Filename of the MAT-file with the results of computations.
% picfile : str
%     Filename of the resulting picture file.
%
r = load(matfile); % Load results from the `matfile` into `r` struct.

figure
set(gcf, 'Position', [200 200 4.5*96 8.3436*96]);
subplot(3, 1, 1, 'align');
hold on
plot(r.grid, r.znd_all.u, '-')
plot(r.grid, r.znd_all.l, '--')
plot(r.grid, r.znd_all.w, '-.')
title('ZND solution')
legend('u', '\lambda', '\omega', 'Location', 'northwest')

subplot(3, 1, 2, 'align');
hold on
plot(r.grid, r.pert(1, :))
plot(r.grid, r.pert(3, :))
title('Real part')
legend('Re u', 'Re \lambda', 'Location', 'southwest');

subplot(3, 1, 3, 'align');
hold on
plot(r.grid, r.pert(2, :))
plot(r.grid, r.pert(4, :))
xlabel('x')
title('Imaginary part')
legend('Im u', 'Im \lambda', 'Location', 'southwest');

export_fig(picfile);
end