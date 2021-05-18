function format_figure( han )
%UNTITLED3 Summary of this function goes here
%   changes the defaults on a figure... like I do for all figures

if nargin == 0
    han = gca;
end

set(han,'linewidth',2,'fontsize',18,'fontweight','bold');

end

