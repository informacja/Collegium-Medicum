function FigArticle(filename, format, fontSize, figSizeScalar)
% FigArticle - Save current figure in Nature-style format
%
% Usage:
%   FigArticle('F1')                    % vector PDF, default font 12
%   FigArticle('F1','image')            % raster image
%   FigArticle('F1','vector',14)        % vector image, font size 14
%
% Inputs:
%   filename - string, name of the file
%   format   - optional, 'vector' (default) or 'image'
%   fontSize - optional, numeric font size (default = 12)

% ---------------- Handle defaults ----------------
if nargin < 2
    format = 'vector';
end
if nargin < 3
    fontSize = 12;  % default font size
end
if nargin < 4
    figSizeScalar = 1;  % default font size
end

inerpreter = "latex";

% ---------------- Ensure Results folder exists ----------------
if ~exist('Results','dir')
    mkdir('Results');
end

% ---------------- Configure current figure ----------------
set(gcf, 'Color', 'w', ...                     % white background
         'Units', 'normalized', ...            % normalized units
         'Position', [1.2 0.1 0.6 0.6]);       % default position/size

% ---------------- Force Helvetica font for axes ----------------
axes_handles = findobj(gcf, 'Type', 'axes');
for ax = axes_handles'
    ax.FontName = 'Helvetica';
    ax.FontSize = fontSize;
    ax.LineWidth = 1.2;
    ax.TickDir = 'out';
    
    % Set labels and title fonts
    ax.XLabel.FontName = 'Helvetica';
    ax.XLabel.FontSize = fontSize;
    ax.XLabel.Interpreter = inerpreter;
    ax.YLabel.FontName = 'Helvetica';
    ax.YLabel.FontSize = fontSize;
    ax.YLabel.Interpreter = inerpreter;
    if(inerpreter=="latex") mod=(strrep(yticklabels,'-','$-$')); ax.YTickLabel = mod; end
    ax.Title.FontName  = 'Helvetica';
    ax.Title.FontSize  = fontSize;    
    ax.Title.Interpreter = inerpreter;
    ax.TickLabelInterpreter = inerpreter;
end

% ---------------- Force Helvetica font for legends ----------------
legend_handles = findobj(gcf, 'Type', 'Legend');
for lg = legend_handles'
    lg.FontName = 'Helvetica';
    lg.FontSize = fontSize - 1;  % slightly smaller for legend
    lg.Interpreter = inerpreter;
end

% ---------------- Build full filename with PDF extension ----------------
filenameFull = fullfile('Results', [filename, '.pdf']);

inches = 7.5; % full Page
width = inches*1.618*figSizeScalar; % golden ratio
hight = inches*figSizeScalar;
set(gcf,'Units','inches');                        
set(gcf,'Position', [0 0 width hight]); 

% ---------------- Export figure ----------------
switch lower(format)
    case 'vector'
        exportgraphics(gcf, filenameFull, 'ContentType','vector', 'BackgroundColor','w');
    case 'image'
        exportgraphics(gcf, filenameFull, 'ContentType','image', 'BackgroundColor','w', 'Resolution',300);
    otherwise
        error('Format must be ''vector'' or ''image''.');
end

fprintf('✅ Figure saved: %s\n', filenameFull);
end
