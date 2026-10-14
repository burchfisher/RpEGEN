function test_rpegen_permutation_gui
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
files = {'permsim.m','rpegen_permutation_profile.m','rpegen_permutation_tab.m', ...
    'rpegen_comparison_tab.m','RpEGEN_GUI.m'};
for file = files
    issues = checkcode(fullfile(root,file{1}),'-id');
    for k = 1:numel(issues), fprintf('%s:%d %s\n',file{1},issues(k).line,issues(k).message); end
end
angles = repelem((0.5:179.5)',12); values = repmat((0:11)',180,1);
a = table(angles,values,'VariableNames',{'Angle','Raw_Data'});
b = table(angles,values+3,'VariableNames',{'Angle','Raw_Data'});
groups = struct('Name',{'Untreated','Treated'},'File',{'test-a.mat','test-b.mat'}, ...
    'Samples',{3,3},'Raw',{a,b});
fig = uifigure('Visible','off','Position',[100 100 1180 780]);
cleanup = onCleanup(@() delete(fig));
grid = uigridlayout(fig,[1 1]); tabs = uitabgroup(grid); tab = uitab(tabs,'Title','Permutation Test');
busyEvents = []; api = rpegen_permutation_tab(tab,root,@trackBusy,@() []);
runButton = findall(fig,'Tag','RunPermutation');
assert(strcmp(runButton.Enable,'off'));
api.SetGroups(groups(1)); assert(strcmp(runButton.Enable,'off'));
api.SetGroups(groups); assert(strcmp(runButton.Enable,'on'));
mode = findall(fig,'Tag','PermutationMode'); mode.Value = 'Manual'; feval(mode.ValueChangedFcn,mode,[]);
draws = findall(fig,'Tag','PermutationIterations'); draws.Value = 257;
width = findall(fig,'Tag','PermutationBin'); width.Value = 10;
parallelToggle = findall(fig,'Tag','PermutationParallel'); parallelToggle.Value = false;
before = rng; api.Run(); [serial,metadata] = api.GetResults();
expected = rpegen_permutation_profile(a,b,10,257);
assert(isequal(serial,expected) && isequal(before,rng));
assert(isequal(busyEvents,[true false]));
assert(strcmp(metadata.RunInfo.Execution,'serial') && api.HasResults());
assert(all(serial.Iterations==257) && all(serial.LegacyP==serial.Exceedances./serial.Iterations));
assert(numel(findall(fig,'Tag','PermutationP'))==1);
assert(strcmp(findall(fig,'Tag','PermutationP').DisplayName,'Permutation p-value'));
steps = findall(fig,'Tag','PermutationP');
assert(isa(steps,'matlab.graphics.chart.primitive.Stair'));
assert(isequal(steps.XData(:),(0:10:180)'));
assert(isequal(steps.YData(:),[serial.LegacyP; serial.LegacyP(end)]));
intervals = findall(fig,'Tag','PermutationIntervals');
assert(~intervals.Value && isempty(findall(fig,'Tag','PermutationMC')));
ax = findall(fig,'Tag','PermutationAxes'); ax.XLim = [20 80]; ax.YLim = [0 0.1];
shading = findall(fig,'Tag','PermutationShading');
assert(shading.Value && numel(findall(fig,'Tag','PermutationRibbon'))==1);
shading.Value = false; feval(shading.ValueChangedFcn,shading,[]);
assert(isempty(findall(fig,'Tag','PermutationRibbon')) && api.HasResults());
assert(isequal(ax.XLim,[20 80]) && isequal(ax.YLim,[0 0.1]));
shading.Value = true; feval(shading.ValueChangedFcn,shading,[]);
assert(numel(findall(fig,'Tag','PermutationRibbon'))==1);
intervals.Value = true; feval(intervals.ValueChangedFcn,intervals,[]);
assert(isequal(ax.XLim,[20 80]) && isequal(ax.YLim,[0 0.1]));
bars = findall(fig,'Tag','PermutationMC');
assert(bars.CapSize>0 && bars.LineWidth>=1.25);
assert(isequal(bars.YNegativeDelta(:),serial.LegacyP-serial.MC95L));
assert(isequal(bars.YPositiveDelta(:),serial.MC95U-serial.LegacyP));
api.ExportResults(fullfile(tempdir,'rpegen_permutation_gui.csv'));
api.ExportResults(fullfile(tempdir,'rpegen_permutation_gui.mat'));
saved = load(fullfile(tempdir,'rpegen_permutation_gui.mat'));
assert(isequal(saved.report,serial) && strcmp(saved.metadata.GroupA,'Untreated'));
csv = readtable(fullfile(tempdir,'rpegen_permutation_gui.csv'));
assert(height(csv)==18 && all(csv.Iterations==257));
focus = findall(fig,'Tag','PermutationFocus'); focus.Value = true; feval(focus.ValueChangedFcn,focus,[]);
ax = findall(fig,'Tag','PermutationAxes'); assert(isequal(ax.YLim,[0 0.12]));
mustFail(@() api.ExportFigure(fullfile(tempdir,'rpegen_permutation_gui.png')),'RpEGEN:FigureFormat');
api.ExportFigure(fullfile(tempdir,'rpegen_permutation_gui.pdf'));
intervals.Value = false; feval(intervals.ValueChangedFcn,intervals,[]);
assert(isempty(findall(fig,'Tag','PermutationMC')) && api.HasResults());
[unchanged,details] = api.GetResults();
assert(isequal(unchanged,serial) && isequal(details,metadata));
assert(strcmp(findall(fig,'Tag','PermutationP').DisplayName,'Permutation p-value'));
api.ExportFigure(fullfile(tempdir,'rpegen_permutation_gui_no_intervals.pdf'));
feval(width.ValueChangedFcn,width,[]); assert(~api.HasResults());
mustFail(@() api.ExportResults(fullfile(tempdir,'must_not_export.csv')),'RpEGEN:Permutation');
groupB = findall(fig,'Tag','PermutationGroupB'); groupB.Value = 'test-a.mat';
feval(groupB.ValueChangedFcn,groupB,[]); assert(strcmp(runButton.Enable,'off'));
groupB.Value = 'test-b.mat'; feval(groupB.ValueChangedFcn,groupB,[]);
api.SetBusy(true); assert(strcmp(runButton.Enable,'off'));
assert(strcmp(intervals.Enable,'off'));
assert(strcmp(shading.Enable,'off'));
api.SetBusy(false); assert(strcmp(runButton.Enable,'on'));
assert(strcmp(intervals.Enable,'on'));
assert(strcmp(shading.Enable,'on'));

% Cancellation during a bin must discard partial results and restore controls.
draws.Value = 1e7;
timerObject = timer('StartDelay',0.2,'TimerFcn',@cancelTest);
timerCleanup = onCleanup(@() delete(timerObject)); start(timerObject);
mustFail(@() api.Run(),'RpEGEN:Cancelled');
assert(~api.HasResults() && strcmp(runButton.Enable,'on'));
assert(contains(findall(fig,'Tag','PermutationStatus').Text,'Cancelled'));
clear timerCleanup;

if license('test','Distrib_Computing_Toolbox')
    draws.Value = 257; parallelToggle.Value = true; api.Run();
    [parallel,metadata] = api.GetResults();
    assert(isequal(parallel,serial) && strcmp(metadata.RunInfo.Execution,'parfeval'));
    legacyParallel = rpegen_permutation_profile(a,b,10,257,struct('UseParallel',true));
    assert(isequal(legacyParallel,serial));
    draws.Value = 1e7;
    timerObject = timer('StartDelay',0.3,'TimerFcn',@cancelTest);
    timerCleanup = onCleanup(@() delete(timerObject)); start(timerObject);
    mustFail(@() api.Run(),'RpEGEN:Cancelled');
    assert(~api.HasResults() && strcmp(runButton.Enable,'on'));
    clear timerCleanup;
end

% Auto defaults, successful rerun after cancellation, and validation failure.
parallelToggle.Value = false; mode.Value = 'Auto'; feval(mode.ValueChangedFcn,mode,[]);
width.Value = 10; api.Run(); [automatic,metadata] = api.GetResults();
assert(all(automatic.Iterations>=10000 & automatic.Iterations<=50000));
assert(metadata.AutoMinIterations==10000 && metadata.AutoMaxIterations==50000);
assert(metadata.AutoHalfWidth==0.002);
assert(metadata.RunInfo.TimeBudgetSeconds==300);
groups(1).Raw = a(a.Angle<10,:); api.SetGroups(struct('Name',{},'File',{},'Samples',{},'Raw',{}));
api.SetGroups(groups); mustFail(@() api.Run(),'RpEGEN:Permutation');
assert(~api.HasResults() && strcmp(runButton.Enable,'on'));
disp('GUI serial/futures/parfor parity, progress, cancellation, rerun, input validation, stale-result clearing, and exports passed.');

    function trackBusy(value)
        busyEvents(end+1) = value;
    end
    function cancelTest(~,~)
        assert(strcmp(runButton.Enable,'off'));
        api.Cancel();
    end
end

function mustFail(callback,identifier)
try, callback(); catch problem, assert(strcmp(problem.identifier,identifier),problem.message); return; end
error('Expected failure did not occur.');
end
