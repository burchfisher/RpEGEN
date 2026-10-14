function verify_rpegen_permutation_integration
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
fig = RpEGEN_GUI; cleanup = onCleanup(@() delete(fig));
compareTab = findall(fig,'Type','uitab','Title','Compare Groups');
compareTab.Parent.SelectedTab = compareTab;
button = findall(fig,'Tag','ExampleGroups'); feval(button.ButtonPushedFcn,button,[]);
views = findall(fig,'Tag','ComparisonViews');
views.SelectedTab = findall(fig,'Type','uitab','Title','Permutation Test');
feval(views.SelectionChangedFcn,views,[]);
assert(strcmp(findall(fig,'Tag','ComparisonGroups').Visible,'off'));
a = findall(fig,'Tag','PermutationGroupA'); b = findall(fig,'Tag','PermutationGroupB');
a.Value = fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat');
b.Value = fullfile(root,'JoVE_dataset','Output','IWR1_4dpi.mat');
feval(a.ValueChangedFcn,a,[]); feval(b.ValueChangedFcn,b,[]);
mode = findall(fig,'Tag','PermutationMode'); mode.Value = 'Manual'; feval(mode.ValueChangedFcn,mode,[]);
draws = findall(fig,'Tag','PermutationIterations'); draws.Value = 1e7;
parallelToggle = findall(fig,'Tag','PermutationParallel'); parallelToggle.Value = false;
runButton = findall(fig,'Tag','RunPermutation'); cancelButton = findall(fig,'Tag','CancelPermutation');
analysisButton = findall(fig,'Type','uibutton','Text','Run Analysis');
busyChecked = false;
t = timer('StartDelay',0.3,'TimerFcn',@cancelTest); timerCleanup = onCleanup(@() delete(t)); start(t);
feval(runButton.ButtonPushedFcn,runButton,[]);
assert(busyChecked && contains(findall(fig,'Tag','PermutationStatus').Text,'Cancelled'));
assert(strcmp(analysisButton.Enable,'on') && strcmp(runButton.Enable,'on'));
clear timerCleanup;

mode.Value = 'Auto'; feval(mode.ValueChangedFcn,mode,[]); parallelToggle.Value = true;
feval(runButton.ButtonPushedFcn,runButton,[]);
results = findall(fig,'Tag','PermutationResults'); report = results.Data;
assert(height(report)==180 && all(report.Iterations>=10000 & report.Iterations<=50000));
assert(all(report.LegacyP==report.Exceedances./report.Iterations));
steps = findall(fig,'Tag','PermutationP');
assert(isequal(steps.XData(:),(0:180)'));
assert(isequal(steps.YData(:),[report.LegacyP; report.LegacyP(end)]));
ribbon = findall(fig,'Tag','PermutationRibbon');
assert(numel(ribbon)==1);
covered = any(report.Angle>=ribbon.XData(1,:) & report.Angle<ribbon.XData(3,:),2);
assert(isequal(covered,report.LegacyP<0.05));
assert(all(ribbon.XData(1,:)>=0 & ribbon.XData(3,:)<=180));
assert(all(ribbon.XData(3,:)>ribbon.XData(1,:)));
status = findall(fig,'Tag','PermutationStatus'); disp(status.Text);
assert(startsWith(status.Text,'Complete:'));
assert(strcmp(analysisButton.Enable,'on') && strcmp(runButton.Enable,'on'));
focus = findall(fig,'Tag','PermutationFocus'); focus.Value = true; feval(focus.ValueChangedFcn,focus,[]);
drawnow; pause(1); drawnow;
exportapp(fig,fullfile(root,'comparison_preview','permutation_gui_example.png'));
writetable(report,fullfile(root,'comparison_preview','permutation_gui_example.csv'));
resultViews = findall(fig,'Tag','PermutationViews');
resultViews.SelectedTab = findall(fig,'Type','uitab','Title','Results Table');
drawnow; pause(0.5); exportapp(fig,fullfile(root,'comparison_preview','permutation_gui_table.png'));
angle = findall(fig,'Tag','PermutationBin'); angle.Value = 5; feval(angle.ValueChangedFcn,angle,[]);
assert(isempty(results.Data));
exportButton = findall(fig,'Type','uibutton','Text','Export Figure...'); assert(strcmp(exportButton.Enable,'off'));
views.SelectedTab = findall(fig,'Type','uitab','Title','Median Profiles');
feval(views.SelectionChangedFcn,views,[]);
assert(strcmp(findall(fig,'Tag','ComparisonGroups').Visible,'on'));
disp('Full GUI: global busy state, real-data Auto run, cancellation, table, plot, view switching, and stale-result protection passed.');

    function cancelTest(~,~)
        busyChecked = strcmp(analysisButton.Enable,'off') && strcmp(runButton.Enable,'off');
        feval(cancelButton.ButtonPushedFcn,cancelButton,[]);
    end
end
