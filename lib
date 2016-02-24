function export_fig_in_pdf(filename, size_inches)
%MY_EXPORT_FIG Wrapper around `export_fig` to tight layout of the image.
% Parameters
% ----------
% filename : string
%     Filename of the image.
% size_inches : array (2x1)
%     Size of the image in inches.
%
size_x = size_inches(1);
size_y = size_inches(2);

% To obtain information about the resolution, see
% http://www.mathworks.com/matlabcentral/answers/100792-in-matlab-how-do-i-obtain-information-about-my-screen-resolution-and-screen-size
res = 96; % Resolution of my screen is 96 dots per inch.

set(gcf, 'Position', [200, 200, size_x * res, size_y * res]);
tight_layout();
export_fig(filename, '-nocrop');
end


%--------------------------------------------------------------------------
function tight_layout
% Tight layout.
% See http://www.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end