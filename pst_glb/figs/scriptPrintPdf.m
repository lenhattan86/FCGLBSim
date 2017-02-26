PS_CMD_FORMAT='ps2pdf -dEmbedAllFonts#true -dSubsetFonts#true -dEPSCrop#false -dPDFSETTINGS#/prepress %s %s';
fig_path = '../pdf_figs/';
%%
epsFiles = {'load_freq';
             'gen_freq'}
%%


%% call ps2pdf
for i=1:length(epsFiles)
    fileName = epsFiles{i};
    epsFile = [fileName '.eps'];
    pdfFile = [fig_path fileName  '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end
%%