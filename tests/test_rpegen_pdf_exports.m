function test_rpegen_pdf_exports
% Exercise publication exports independently of extraction and long simulations.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
output = fullfile(tempdir,'rpegen_pdf_exports');
if ~isfolder(output), mkdir(output); end
fig = figure('Visible','off','Color','w'); cleanup = onCleanup(@() close(fig));
plot(1:5,[1 3 2 5 4]); xlabel('Editable angle'); ylabel('Editable width');
file = rpegen_export_figure(fig,fullfile(output,'extension_added'));
assert(strcmp(file,fullfile(output,'extension_added.pdf')) && isfile(file));
rpegen_export_figure(fig,fullfile(output,'uppercase.PDF'));
try
    rpegen_export_figure(fig,fullfile(output,'rejected.png'));
    error('Test:ExpectedError','PNG was unexpectedly accepted.');
catch problem
    assert(strcmp(problem.identifier,'RpEGEN:FigureFormat'));
end
assert(~isfile(fullfile(output,'rejected.png')));
clear cleanup;

fig = uifigure('Visible','off','Position',[100 100 1180 780]);
cleanup = onCleanup(@() delete(fig));
grid = uigridlayout(fig,[1 1]); tabs = uitabgroup(grid);
comparisonTab = uitab(tabs,'Title','Compare Groups');
comparison = rpegen_comparison_tab(comparisonTab,root);
source = fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat');
comparison.AddFiles({source});
comparison.Export(fullfile(output,'heatmaps.pdf'),'Heatmaps');
comparison.Export(fullfile(output,'median_profiles.pdf'),'Median Profiles');

angles = repelem((5:10:175)',12); values = repmat((0:11)',18,1);
a = table(angles,values,'VariableNames',{'Angle','Raw_Data'});
b = table(angles,values+3,'VariableNames',{'Angle','Raw_Data'});
groups = struct('Name',{'Untreated','Treated'},'File',{'a.mat','b.mat'}, ...
    'Samples',{3,3},'Raw',{a,b});
comparison.Permutation.SetGroups(groups);
mode = findall(fig,'Tag','PermutationMode'); mode.Value = 'Manual';
feval(mode.ValueChangedFcn,mode,[]);
iterations = findall(fig,'Tag','PermutationIterations'); iterations.Value = 200;
bin = findall(fig,'Tag','PermutationBin'); bin.Value = 10;
parallel = findall(fig,'Tag','PermutationParallel'); parallel.Value = false;
comparison.Permutation.Run();
comparison.Export(fullfile(output,'permutation.pdf'),'Permutation Test');

for name = {'heatmaps','median_profiles','permutation'}
    file = fullfile(output,[name{1} '.pdf']);
    fid = fopen(file,'r'); fileCleanup = onCleanup(@() fclose(fid));
    assert(strcmp(fread(fid,5,'*char')','%PDF-'));
    clear fileCleanup;
end
fprintf('PDF export checks passed. Inspect vector content and rendering in %s\n',output);
end
