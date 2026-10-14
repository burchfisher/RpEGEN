function test_rpegen_confidence
% Verify legacy saved limits and their display, without changing result files.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
files = fullfile(root,'JoVE_dataset','Output', ...
    {'DMSO_4dpi.mat','IWR1_4dpi.mat','DMSO_9dpf.mat','IWR1_9dpf.mat'});
binCount = 0;
for f = 1:numel(files)
    loaded = load(files{f}); name = fieldnames(loaded); data = loaded.(name{1});
    for width = [1 5 10]
        bins = data.(sprintf('BIN_%d_deg_all',width));
        for b = 1:height(bins)
            rd = sort(data.RAW_data.Raw_Data(data.RAW_data.Angle >= (b-1)*width & ...
                data.RAW_data.Angle < b*width));
            n = numel(rd);
            j = round(n/2-1.96*sqrt(n/4)); k = round(n/2+1.96*sqrt(n/4));
            assert(j >= 1 && k <= n,'Example data contains a sparse bin.');
            assert(isequal([double(rd(j)),double(rd(k)),double(median(rd))], ...
                [bins.CI95L(b),bins.CI95U(b),bins.Median(b)]));
            binCount = binCount+1;
        end
    end
end
fprintf('Reproduced saved legacy medians and limits in %d bins.\n',binCount);

fig = uifigure('Visible','off','Position',[100 100 1180 780]);
cleanup = onCleanup(@() delete(fig));
tabs = uitabgroup(fig,'Position',[10 10 1160 760]);
tab = uitab(tabs,'Title','Compare Groups');
api = rpegen_comparison_tab(tab,root); api.AddFiles(files(1:2));
angle = findall(fig,'Tag','HeatmapAngleBin');
toggle = findall(fig,'Tag','MedianConfidenceIntervals');
bootstrapToggle = findall(fig,'Tag','BootstrapConfidenceIntervals');
assert(toggle.Value && bootstrapToggle.Value);
assert(strcmp(toggle.Text,'Pixel-based median CI (nominal 95%)'));
assert(~contains(lower(toggle.Tooltip),'legacy'));
ax = findall(fig,'Tag','GroupMedianAxes');
assert(isscalar(ax) && isempty(findall(fig,'Tag','GroupBootstrapAxes')));
for width = [1 5 10 5 1]
    angle.Value = width; feval(angle.ValueChangedFcn,angle,[]);
    assert(isequal(ax,findall(fig,'Tag','GroupMedianAxes')));
    bootstrapBands = findall(ax,'Tag','BootstrapConfidenceBand');
    assert(numel(bootstrapBands)==2);
    bands = findall(ax,'Tag','MedianConfidenceBand'); assert(numel(bands)==2);
    assert(all([bands.FaceAlpha]>max([bootstrapBands.FaceAlpha])));
    curves = findall(ax,'Type','line'); assert(numel(curves)==2);
    for k = 1:2
        loaded = load(files{k}); name = fieldnames(loaded);
        bins = loaded.(name{1}).(sprintf('BIN_%d_deg_all',width));
        band = bands(arrayfun(@(b) strcmp(b.UserData,name{1}),bands));
        assert(isequal(band.XData(:),[bins.Angle;flipud(bins.Angle)]));
        assert(isequal(band.YData(:),[bins.CI95L;flipud(bins.CI95U)]));
        curve = curves(arrayfun(@(c) startsWith(c.DisplayName,name{1}),curves));
        assert(isequaln(curve.YData(:),bins.Median));
    end
end
beforeToggle = rng;
xlim(ax,[40 60]); ylim(ax,[50 170]);
for state = [1 0;0 1;0 0;1 1]'
    toggle.Value = logical(state(1)); feval(toggle.ValueChangedFcn,toggle,[]);
    bootstrapToggle.Value = logical(state(2)); feval(bootstrapToggle.ValueChangedFcn,bootstrapToggle,[]);
    assert(numel(findall(ax,'Tag','MedianConfidenceBand'))==2*state(1));
    assert(numel(findall(ax,'Tag','BootstrapConfidenceBand'))==2*state(2));
    assert(numel(findall(ax,'Type','line'))==2 && isequal(ax.XLim,[40 60]) && isequal(ax.YLim,[50 170]));
    assert(~contains(lower(findall(fig,'Tag','MedianCIStatus').Text),'legacy'));
    api.Export(fullfile(tempdir,sprintf('rpegen_ci_%d_%d.pdf',state)),'Median Profiles');
end
assert(isequal(beforeToggle,rng),'CI toggles or exports changed the RNG.');
groupTable = findall(fig,'Tag','ComparisonGroups');
for k = 1:2
    feval(groupTable.CellEditCallback,groupTable,struct('Indices',[k 1],'NewData',false));
end
ax = findall(fig,'Tag','GroupMedianAxes');
assert(isempty(findall(ax,'Tag','MedianConfidenceBand')));
assert(isempty(findall(ax,'Tag','BootstrapConfidenceBand')));
assert(isempty(findall(ax,'Type','line')));
for k = 1:2
    feval(groupTable.CellEditCallback,groupTable,struct('Indices',[k 1],'NewData',true));
end
api.SetBusy(true); assert(strcmp(toggle.Enable,'off') && strcmp(bootstrapToggle.Enable,'off'));
api.SetBusy(false); assert(strcmp(toggle.Enable,'on') && strcmp(bootstrapToggle.Enable,'on'));
ax = findall(fig,'Tag','GroupMedianAxes');
assert(numel(findall(ax,'Tag','MedianConfidenceBand'))==2);
angle.Value = 5; feval(angle.ValueChangedFcn,angle,[]);
api.AddFiles(files(3:4));
assert(numel(findall(ax,'Type','line'))==4);
assert(numel(findall(ax,'Tag','MedianConfidenceBand'))==4);
assert(numel(findall(ax,'Tag','BootstrapConfidenceBand'))==4);
views = findall(fig,'Tag','ComparisonViews');
views.SelectedTab = findall(fig,'Type','uitab','Title','Median Profiles');
fig.Visible = 'on'; drawnow; pause(1); drawnow;
exportapp(fig,fullfile(root,'comparison_preview','overlaid_ci_gui.png'));
api.Export(fullfile(root,'comparison_preview','overlaid_ci_export.pdf'),'Median Profiles');
api.Export(fullfile(tempdir,'rpegen_overlaid_ci.pdf'),'Median Profiles');
disp('One axes, unchanged medians, both independent CI toggles, zoom preservation, selection, busy state, and exports passed.');
end
