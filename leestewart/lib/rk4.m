function [sol] = rk4(fh, times, ic)
%RK4 Summary of this function goes here
%   Detailed explanation goes here
% y must be length(ic) \times length(times).

assert(isvector(ic))

[~, n] = size(ic);

if n == 1  % column vector
    y = ic;
else       % row vector, must be transposed
    y = ic';
end
    
sol = zeros(length(y), length(times));
sol(:, 1) = y;
    
for i = 2:length(times)
    curr = i-1;  % index of time layer on which we know data.
    t = times(curr);
    h = times(i) - t;
    h2 = 0.5 * h;
    
    k1 = fh(t, y, curr);
    k2 = fh(t + h2, y + h2*k1, curr);
    k3 = fh(t + h2, y + h2*k2, curr);
    k4 = fh(t + h,  y + h *k3, curr);
    y = y + (h/6) * (k1 + 2*(k2 + k3) + k4);
    sol(:, i) = y;
end
end