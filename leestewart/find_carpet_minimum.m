function [xmin, ymin] = find_carpet_minimum(H)
% Find minimum of the function `H`.
% 
% As function `H` is complex-valued, the minimum is found in an absolute
% value. Indices of the row and column of the minimum are returned.
% 
% Parameters
% ----------
% H : matrix
%     Complex-valued function `H` stemming from a carpet search.
% 
% Returns
% -------
% xmin, ymin : int
%     Indices of the row and the column, resp., of the minimum.
[M, I] = min(abs(H));

[mm, ii] = min(M);

xmin = I(ii);
ymin = ii;

end
