function verify_rpegen_update
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
fig = RpEGEN_GUI('Visible','on'); cleanup = onCleanup(@() delete(fig));
tab = findall(fig,'Type','uitab','Title','Compare Groups'); tab.Parent.SelectedTab = tab;
button = findall(fig,'Tag','ExampleGroups'); feval(button.ButtonPushedFcn,button,[]);
views = findall(fig,'Tag','ComparisonViews');
views.SelectedTab = findall(fig,'Type','uitab','Title','Median Profiles');
for tag = {'MedianConfidenceIntervals','BootstrapConfidenceIntervals'}
    toggle = findall(fig,'Tag',tag{1});
    toggle.Value = false; feval(toggle.ValueChangedFcn,toggle,[]);
    toggle.Value = true; feval(toggle.ValueChangedFcn,toggle,[]);
end
drawnow; pause(1); drawnow;
exportapp(fig,fullfile(root,'comparison_preview','main_gui_ci.png'));
pause(1); drawnow;
exportapp(fig,fullfile(root,'comparison_preview','main_gui_ci.png'));
delete(fig); clear cleanup;
a = load(fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'));
b = load(fullfile(root,'JoVE_dataset','Output','IWR1_4dpi.mat'));
[report,runInfo] = rpegen_permutation_profile(a.DMSO_4dpi.RAW_data, ...
    b.IWR1_4dpi.RAW_data,1,'auto',struct('UseParallel',true,'TimeBudgetSeconds',300));
assert(height(report)==180 && all(report.Iterations>=10000 & report.Iterations<=50000));
save(fullfile(root,'comparison_preview','permutation_auto_example.mat'),'report','runInfo');
writetable(report,fullfile(root,'comparison_preview','permutation_auto_example.csv'));
disp(runInfo);
fprintf('Auto example: %d-%d draws per bin, %d MC intervals crossing 0.05.\n', ...
    min(report.Iterations),max(report.Iterations),nnz(report.NearThreshold));
end
