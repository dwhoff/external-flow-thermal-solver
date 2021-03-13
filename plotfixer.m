% Plotfixer

plotlsize = 2; %thickness of plotted lines, in points
axislsize = 1.5; %thickness of tick marks and borders, in points
markersize = 8;  %size of line markers, default is 6

axisfont = 'Helvetica'; % changes appearance of axis numbers
axisfontsize = 18;            % in points
axisfontweight = 'normal';    % options are 'light' 'normal' 'demi' 'bold' 
axisfontitalics = 'normal';   % options are 'normal' 'italic' 'oblique'

legendfontsize = 20; % changes text in the legend

labelfont = 'Helvetica';  % changes x, y, and z axis labels
labelfontsize = 24;  
labelfontweight = 'normal'; 
labelfontitalics = 'normal';

titlefont = 'Helvetica';  % changes title
titlefontsize = 18;
titlefontitalics = 'normal';

textfont = 'Helvetica';   % changes text
textfontsize = 18;
textfontweight = 'normal';


%stop changing things below this line
%----------------------------------------------------
axesh = findobj('Type', 'axes');
lineh = findobj(axesh, 'Type', 'line');
axestexth = findobj(axesh, 'Type', 'text');

set(lineh, 'LineWidth', plotlsize)
set(lineh, 'MarkerSize', markersize)
set(axesh, 'LineWidth', axislsize)
set(axesh, 'FontName', axisfont)
set(axesh, 'FontSize', axisfontsize)
set(axesh, 'FontWeight', axisfontweight)
set(axesh, 'FontAngle', axisfontitalics)
set(axestexth, 'FontName', textfont)
set(axestexth, 'FontSize', textfontsize)
set(axestexth, 'FontWeight', textfontweight)
set(axesh, 'Box','on')
for i = 1:1:size(axesh)
    axes(axesh(i))
    set(get(gca,'XLabel'), 'FontName', labelfont)
    set(get(gca,'XLabel'), 'FontSize', labelfontsize)
    set(get(gca,'XLabel'), 'FontWeight', labelfontweight)
    set(get(gca,'XLabel'), 'FontAngle', labelfontitalics)
    set(get(gca,'YLabel'), 'FontName', labelfont)
    set(get(gca,'YLabel'), 'FontSize', labelfontsize)
    set(get(gca,'YLabel'), 'FontWeight', labelfontweight)
    set(get(gca,'YLabel'), 'FontAngle', labelfontitalics)
    set(get(gca,'ZLabel'), 'FontName', labelfont)
    set(get(gca,'ZLabel'), 'FontSize', labelfontsize)
    set(get(gca,'ZLabel'), 'FontWeight', labelfontweight)
    set(get(gca,'ZLabel'), 'FontAngle', labelfontitalics)
    set(get(gca,'Title'), 'FontName', titlefont)
    set(get(gca,'Title'), 'FontSize', titlefontsize)
    set(get(gca,'Title'), 'FontAngle', titlefontitalics)
    set(get(gca,'Legend'), 'FontSize', legendfontsize)
end

clear axesh axestexth axisfont axisfontitalics axisfontsize axisfontweight ...
    axislsize i labelfont labelfontitalics labelfontsize labelfontweight ...
    legendfont legendfontitalics legendfontsize legendfontweight legendh ...
    lineh markersize plotlsize textfont textfontitalics textfontsize ...
    textfontweight titlefont titlefontitalics titlefontsize titlefontweight