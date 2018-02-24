%% path

OUTPUT_FOLDER = 'output/';
LOCAL_FIG = 'figs/';
%fig_path = 'C:\Users\NhatTan\Dropbox\Papers\FCGLB\FCGLB\cdc_17_paper\images\';
fig_path = 'figs/';
PS_CMD_FORMAT='ps2pdf -dEmbedAllFonts#true -dSubsetFonts#true -dEPSCrop#false -dPDFSETTINGS#/prepress %s %s';
%% String
strProposed = 'proposed';
strOLC = 'OLC';
strOptimal = 'Lower bound';
strSubOptimal = 'suboptimal';
strFCOptimal = 'opt';

strFCLocal = 'ind';

strNone = 'none';

strDelay = 'delay (sec)';

strFrequency = 'frequency (Hz)';
strTime = 'time (seconds)';
strCost = 'cost';

strTotalCost = 'total cost (USD/sec)';

strGlobalCost = 'interdependent';
strLocalCost = 'independent';

strPower = 'Power (MW)';

strGamma = '\gamma';

%%
fontSize=10;

fontAxis = fontSize;
fontTitle = fontSize;
fontLegend = fontSize;

FontSize = fontSize;
is_printed = true;
axisWidth = 1.5;

fontName = 'Arial';

%%

figOneCol = [0 0 5 3];
figHalfCol = 1/2*figOneCol;
figTwoThirdCol = 2/3*figOneCol;

%% lines
lineWidth = 1.5;

lineProposed = '-';
lineOLC = '-.';
lineOptimal = '--';
lineNone = ':';

%% bars
barWidth = 0.7;

%% colors
%http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
colorProposed = [237    125    49]/255;
colorOLC = [165    165    165]/255;
colorNone = [68    71   196]/255;
colorOptimal = [0    0   0]/255;

colorLocal = [hex2dec('2c')    hex2dec('7f')    hex2dec('b8')]/255;
colorGlobal = [hex2dec('7f')    hex2dec('cd')    hex2dec('bb')]/255;

