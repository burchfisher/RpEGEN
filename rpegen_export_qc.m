function files = rpegen_export_qc(data,folder,progress)
% Regenerate QC PDFs from saved data without extraction or MAT-file changes.
if nargin<3, progress = @(~,~) []; end
assert(isfolder(folder),'RpEGEN:Output','QC output folder does not exist.');
fig = figure('Visible','off','Color','w','Position',[100 100 1400 850], ...
    'Tag','RpEGENExport');
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
axesHandles = gobjects(1,3);
for j = 1:3, axesHandles(j) = nexttile(layout,j); end
profile = nexttile(layout,4,[1 3]);
files = cell(numel(data.IMG),1);
for k = 1:numel(files)
    progress((k-1)/numel(files),sprintf('Exporting QC PDF %d / %d',k,numel(files)));
    rpegen_draw_qc(axesHandles,profile,data,k);
    files{k} = fullfile(folder,[data.IMGname{k} '_QC_plots.pdf']);
    rpegen_export_figure(fig,files{k});
end
progress(1,sprintf('Exported %d QC PDFs',numel(files)));
end
