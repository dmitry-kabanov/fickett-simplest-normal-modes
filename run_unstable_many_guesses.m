% Test find_roots with steady-piston rear boundary condition.
% Case is unstable ZND solution.
close all; clear all; clc

resultdir = 'results/2016-03-27-towards-shock-many-guesses/';
filename = strcat(resultdir, 'unstable');

% Physical free parameters.
q = 1.7; theta = 2.4;

% Grid resolution;
N = 1e4;

N_EXP = 5;

u_re = linspace(-1, 1, N_EXP);
u_im = linspace(-1, 1, N_EXP);

for k = 1:N_EXP
    filename_k = sprintf('%s_%02d', filename, k);
    % Guess for eigenvalue.
    guess = [u_re(k); u_im(k); 0; 0; 0.0327; 0.6626];

    [params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess);
    save(strcat(filename_k, '.mat'));

    % Plotting part.
    plot_znd_and_perturbations(filename_k);
end

for k = 1:N_EXP
    filename_k = sprintf('%s_%02d', filename, k);
    data(k) = load(strcat(filename_k, '.mat'));    
end

figure
hold on
for k = 1:N_EXP
    plot(data(k).grid, data(k).pert(1, :));
    msg = sprintf('alpha_re = %.9f, alpha_im = %.9f\n', ...
                  data(k).result.root(5), data(k).result.root(6));
    disp(msg);
end
set(gcf, 'Position', [200 200 6*96 3.7082*96]);
xlabel('x');
ylabel('Real part of u');
tight_layout
export_fig(strcat(filename, '_u_re_all.pdf'));

figure
hold on
for k = 1:N_EXP
    plot(data(k).grid, data(k).pert(2, :));
end
set(gcf, 'Position', [200 200 6*96 3.7082*96]);
xlabel('x');
ylabel('Imaginary part of u');
tight_layout
export_fig(strcat(filename, '_u_im_all.pdf'));