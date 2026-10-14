function test_rpegen_qc_navigation
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
fig = RpEGEN_GUI('Visible','off');
cleanup = onCleanup(@() delete(fig));
previous = findall(fig,'Tag','PreviousSample');
next = findall(fig,'Tag','NextSample');
assert(strcmp(previous.Enable,'off') && strcmp(next.Enable,'off'));
clear cleanup;

fig = RpEGEN_GUI('Visible','off','ResultsFile', ...
    fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'));
cleanup = onCleanup(@() delete(fig));
previous = findall(fig,'Tag','PreviousSample');
next = findall(fig,'Tag','NextSample');
menu = findall(fig,'Tag','SampleMenu');
counter = findall(fig,'Tag','SampleCounter');
if ismac
    assert(isscalar(findall(fig,'Text','Show QC PDF in Finder')));
else
    assert(isscalar(findall(fig,'Text','Open QC PDF')));
end
drawnow;
assert(previous.Parent==next.Parent);
assert(previous.Layout.Column==1 && next.Layout.Column==2);
assert(menu.Layout.Column==3 && counter.Layout.Column==4);
assert(strcmp(previous.Enable,'off') && strcmp(next.Enable,'on'));
assert(strcmp(counter.Text,'1 / 11'));
feval(next.ButtonPushedFcn,next,[]);
assert(strcmp(menu.Value,menu.Items{2}) && strcmp(counter.Text,'2 / 11'));
assert(strcmp(previous.Enable,'on'));
feval(previous.ButtonPushedFcn,previous,[]);
assert(strcmp(menu.Value,menu.Items{1}) && strcmp(previous.Enable,'off'));
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','rightarrow'));
assert(strcmp(menu.Value,menu.Items{2}));
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','leftarrow'));
assert(strcmp(menu.Value,menu.Items{1}));
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','leftarrow'));
assert(strcmp(menu.Value,menu.Items{1}));
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{'shift'}},'Key','rightarrow'));
assert(strcmp(menu.Value,menu.Items{1}));
fig.CurrentObject = menu;
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','rightarrow'));
assert(strcmp(menu.Value,menu.Items{1}));
fig.CurrentObject = next;
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','pagedown'));
assert(strcmp(menu.Value,menu.Items{2}));
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','pageup'));
assert(strcmp(menu.Value,menu.Items{1}));
qc = findall(fig,'Type','uitab','Title','Quality Control');
qc.Parent.SelectedTab = findall(fig,'Type','uitab','Title','Results');
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','rightarrow'));
assert(strcmp(menu.Value,menu.Items{1}));
qc.Parent.SelectedTab = qc;
menu.Value = menu.Items{end};
feval(menu.ValueChangedFcn,menu,[]);
assert(strcmp(counter.Text,'11 / 11') && strcmp(next.Enable,'off'));
feval(fig.WindowKeyPressFcn,fig,struct('Modifier',{{}},'Key','rightarrow'));
assert(strcmp(counter.Text,'11 / 11'));
feval(previous.ButtonPushedFcn,previous,[]);
assert(strcmp(counter.Text,'10 / 11') && strcmp(next.Enable,'on'));
menu.Value = menu.Items{1};
feval(menu.ValueChangedFcn,menu,[]);
fig.Visible = 'on'; drawnow; pause(1); drawnow;
exportapp(fig,fullfile(tempdir,'rpegen_qc_navigation.png'));
disp('QC navigation passed: adjacent arrows, dropdown, counter, boundaries, and keyboard.');
end
