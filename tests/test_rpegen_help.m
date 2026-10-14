function test_rpegen_help
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
topics = {'Setup','Quality Control','Results','Heatmaps','Median Profiles','Permutation Test'};
for k = 1:numel(topics)
    [heading,body] = rpegen_help(topics{k});
    assert(contains(heading,topics{k}) && iscell(body) && numel(body)>5);
end
[~,body] = rpegen_help('Results');
for phrase = {'mean 40 but median 20','STD / sqrt(Numel)','MEAN is estimated','MEDIAN','pixels from one eye are related'}
    assert(any(contains(body,phrase{1})));
end
[~,body] = rpegen_help('Heatmaps');
assert(any(contains(body,'200 of 1,000 pixels')) && any(contains(body,'Permutation Test for p-values')));
fig = RpEGEN_GUI('Visible','on'); cleanup = onCleanup(@() delete(fig));
assert(isempty(findall(fig,'Type','uitab','Title','Widths')));
setup = findall(fig,'Type','uitab','Title','Setup'); main = setup.Parent;
assert(isequal(sort(string({main.Children.Title})),sort(["Setup","Quality Control","Results","Compare Groups"])));
button = findall(fig,'Tag','ContextHelp');
for title = {'Setup','Quality Control','Results'}
    main.SelectedTab = findall(fig,'Type','uitab','Title',title{1});
    feval(button.ButtonPushedFcn,button,[]);
    help = findall(groot,'Tag','RpEGENHelp'); assert(isscalar(help) && contains(help.Name,title{1}));
end
example = findall(fig,'Tag','ExampleGroups'); feval(example.ButtonPushedFcn,example,[]);
qcMenu = findall(fig,'Tag','SampleMenu'); assert(strcmp(qcMenu.Enable,'off'));
main.SelectedTab = findall(fig,'Type','uitab','Title','Compare Groups');
views = findall(fig,'Tag','ComparisonViews');
for title = {'Heatmaps','Median Profiles','Permutation Test'}
    views.SelectedTab = findall(fig,'Type','uitab','Title',title{1});
    feval(button.ButtonPushedFcn,button,[]);
    assert(contains(help.Name,title{1}));
    if strcmp(title{1},'Median Profiles')
        control = findall(help,'Tag','RpEGENHelpText');
        [~,expected] = rpegen_help('Median Profiles');
        assert(isequal(control.Value,expected));
        for phrase = {'SAME median','1,000 independent observations','10 independent eyes', ...
                'A, A, C','original angles','Pixels within an eye are related', ...
                'median and permutation results stay the same'}
            assert(any(contains(control.Value,phrase{1})));
        end
    elseif strcmp(title{1},'Permutation Test')
        control = findall(help,'Tag','RpEGENHelpText');
        [~,expected] = rpegen_help('Permutation Test');
        assert(isequal(control.Value,expected));
        for phrase = {'45,618','480 / 10,000 = 0.048','0.044-0.052', ...
                '0.048-0.052','0.05 is still the significance threshold', ...
                'minimum of 10,000','ceiling, not the target','before calculating the final p-value', ...
                'does not keep checking','More pixels','sampling different eyes','fixed seed', ...
                'positive upper limit','false-positive','One random shuffle'}
            assert(any(contains(control.Value,phrase{1})));
        end
        [exampleP,limits] = binofit(480,10000,.05);
        assert(exampleP==.048 && isequal(round(limits,3),[.044 .052]));
        plan = rpegen_permutation_plan;
        assert(plan.Iterations==45618 && plan.ApproximateHalfWidth95==.002);
    end
end
text = findall(help,'Tag','RpEGENHelpText'); assert(any(contains(text.Value,'45,618')));
fig.Visible = 'on'; drawnow; pause(.5); drawnow;
exportapp(fig,fullfile(tempdir,'rpegen_intensity_scope.png'));
views.SelectedTab = findall(fig,'Type','uitab','Title','Median Profiles');
feval(button.ButtonPushedFcn,button,[]);
help.Visible = 'on'; figure(help); drawnow; pause(3); drawnow;
exportapp(help,fullfile(tempdir,'rpegen_ci_help.png'));
views.SelectedTab = findall(fig,'Type','uitab','Title','Permutation Test');
feval(button.ButtonPushedFcn,button,[]); drawnow; pause(3); drawnow;
exportapp(help,fullfile(tempdir,'rpegen_simulation_help.png'));
delete(fig); assert(~isvalid(help)); clear cleanup;
disp('RpEGEN intensity-only scope, all contextual help views, example loading, and help cleanup passed.');
end
