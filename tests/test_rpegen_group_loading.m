function test_rpegen_group_loading
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
fig = uifigure('Visible','on','Position',[100 100 1180 780]); cleanup = onCleanup(@() delete(fig));
grid = uigridlayout(fig,[1 1]); tabs = uitabgroup(grid); tab = uitab(tabs,'Title','Compare Groups');
events = []; captured = false;
api = rpegen_comparison_tab(tab,root,@busyChanged);
source = fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'); api.SetCurrent(source);
example = findall(fig,'Tag','ExampleGroups');
assert(strcmp(example.Interruptible,'off') && strcmp(example.BusyAction,'cancel'));
api.AddFiles({source,source});
assert(numel(api.GetFiles())==1 && isequal(events,[true false]));
assert(numel(findall(fig,'Tag','GroupHeatmap'))==1);
api.AddFiles({source}); assert(isequal(events,[true false]));
feval(example.ButtonPushedFcn,example,[]);
assert(numel(api.GetFiles())==4 && isequal(events,[true false true false]));
assert(strcmp(example.Text,'Example Groups') && captured);
axesBefore = findall(fig,'Tag','GroupHeatmap');
feval(example.ButtonPushedFcn,example,[]);
assert(isequal(axesBefore,findall(fig,'Tag','GroupHeatmap')));
assert(isequal(events,[true false true false]),'Duplicate clicks repeated plotting.');
before = api.GetFiles();
invalid = fullfile(tempdir,'rpegen_invalid_group_loading.mat'); unrelated = 1; save(invalid,'unrelated');
try, api.AddFiles({invalid}); error('Test:ExpectedError','Invalid group accepted.');
catch problem, assert(strcmp(problem.identifier,'RpEGEN:Comparison')); end
assert(isequal(before,api.GetFiles()) && isequal(events,[true false true false true false]));
assert(strcmp(findall(fig,'Tag','ComparisonGroups').Enable,'on'));
for tag = {'AddComparisonGroups','AddCurrentComparisonGroup','ExampleGroups','HeatmapAngleBin'}
    assert(strcmp(findall(fig,'Tag',tag{1}).Enable,'on'));
end
assert(~isempty(findall(fig,'Tag','BootstrapConfidenceBand')));
fprintf('Group loading passed: immediate busy state, repeat/reentrant click guards, deduplication, four example groups, confidence bands, and error recovery.\n');

    function busyChanged(value)
        events(end+1) = value;
        if ~value, return; end
        for tag = {'AddComparisonGroups','AddCurrentComparisonGroup','ExampleGroups','HeatmapAngleBin'}
            assert(strcmp(findall(fig,'Tag',tag{1}).Enable,'off'));
        end
        assert(strcmp(findall(fig,'Tag','ComparisonGroups').Enable,'off'));
        existing = api.GetFiles();
        api.AddFiles({source}); feval(example.ButtonPushedFcn,example,[]);
        assert(isequal(existing,api.GetFiles()));
        if ~captured && strcmp(example.Text,'Loading...')
            drawnow nocallbacks;
            exportapp(fig,fullfile(tempdir,'rpegen_loading_groups.png'));
            captured = true;
        end
    end
end
