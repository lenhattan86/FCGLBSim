%% path

OUTPUT_FOLDER = 'output/';
LOCAL_FIG = 'figs/';
% fig_path = '../../FCGLB/full_paper_5/images/';
fig_path = '../../FCGLB/cdc_17_paper/images/';

PS_CMD_FORMAT='ps2pdf -dEmbedAllFonts#true -dSubsetFonts#true -dEPSCrop#false -dPDFSETTINGS#/prepress %s %s';
%% String
strProposed = 'proposed';
strOLC = 'OLC';
strOptimal = 'Optimal';
strSubOptimal = 'suboptimal';
strFCOptimal = 'FC optimal';
strFCLocal = 'FC independent';

strNone = 'none';

strFrequency = 'frequency (Hz)';
strTime = 'time (seconds)';
strCost = 'cost';

strGlobalCost = 'interdependent';
strLocalCost = 'independent';

%%
fontSize=12;

fontAxis = fontSize;
fontTitle = fontSize;
fontLegend = fontSize;

FontSize = fontSize;
is_printed = true;
axisWidth = 1.5;

fontName = 'Arial';

%%

figOneCol = [0 0 5 3];

%% lines
lineWidth = 1.5;

lineProposed = '-';
lineOLC = '-.';
lineOptimal = '--';
lineNone = ':';

%% bars
barWidth = 0.7;

%% colors
colorProposed = [237    125    49]/255;
colorOLC = [165    165    165]/255;
colorNone = [68    71   196]/255;
colorOptimal = [0    0   0]/255;

