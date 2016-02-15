% Unit test for rk4.m.
% Check that problem of simple (harmonic) oscillator can be integrated
% with several digits of accuracy.
close all
clear all
clc

k = 4.0;
m = 1.0;
omega = sqrt(k / m);

y0 = 1.0;
ydot0 = 0.1;

f = @(t, y, idx) [y(2); -omega^2 * y(1)];

times = linspace(0, 1000*pi, 50000);
sol = rk4(f, times, [y0 ydot0]);