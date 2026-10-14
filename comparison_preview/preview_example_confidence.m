function fig = preview_example_confidence(varargin)
% Standalone comparison preview; original GUI, measurements, and MATs unchanged.
options = inputParser;
addParameter(options,'Visible','on'); parse(options,varargin{:});
root = fileparts(fileparts(mfilename('fullpath')));
names = {'DMSO_9dpf','IWR1_9dpf','DMSO_4dpi','IWR1_4dpi'};
groups = struct('Name',{},'Eyes',{},'Angle',{},'Median',{},'Legacy',{},'Bootstrap',{});
iterations = 5000;
for g = 1:numel(names)
    loaded = load(fullfile(root,'JoVE_dataset','Output',[names{g} '.mat']),names{g});
    data = loaded.(names{g});
    fprintf('Calculating whole-eye bootstrap for %s...\n',names{g});
    ci = bootstrapBounds(data,iterations,20261012+g);
    bins = data.BIN_5_deg_all;
    groups(g) = struct('Name',names{g},'Eyes',numel(data.IMGname), ...
        'Angle',bins.Angle,'Median',bins.Median, ...
        'Legacy',[bins.CI95L bins.CI95U],'Bootstrap',ci);
    fprintf('  %d eyes; median interval width: legacy %.1f, bootstrap %.1f\n', ...
        groups(g).Eyes,median(diff(groups(g).Legacy,1,2)),median(diff(ci,1,2)));
end
colors = lines(numel(groups)); included = true(1,numel(groups));
fig = uifigure('Name','RpEGEN - CI comparison preview','Position',[100 100 1480 850], ...
    'Visible',options.Results.Visible,'Tag','RpEGENConfidencePreview');
layout = uigridlayout(fig,[5 1],'RowHeight',{30,130,34,'1x',40}, ...
    'Padding',[14 12 14 12]);
uilabel(layout,'Text','Compare Groups / Median Profiles / Preview','FontSize',18,'FontWeight','bold');
rows = cell(numel(groups),3);
for g = 1:numel(groups), rows(g,:) = {true,groups(g).Name,groups(g).Eyes}; end
groupTable = uitable(layout,'Data',rows,'ColumnName',{'Show','Group','Animals / eyes'}, ...
    'ColumnEditable',[true false false],'ColumnFormat',{'logical','char','numeric'}, ...
    'ColumnWidth',{60,260,120},'CellEditCallback',@editGroup,'Tag','PreviewGroups');
toolbar = uigridlayout(layout,[1 5],'ColumnWidth',{100,280,150,'1x',140},'Padding',[0 0 0 0]);
uilabel(toolbar,'Text','Confidence bands');
methodMenu = uidropdown(toolbar,'Items',{'Both (side by side)','Pixel-based median CI','Whole-eye bootstrap CI'}, ...
    'ItemsData',{'both','legacy','bootstrap'},'Value','both', ...
    'ValueChangedFcn',@refresh,'Tag','PreviewCIMethod');
uilabel(toolbar,'Text','5-degree bins');
uilabel(toolbar,'Text','Pooled-pixel medians');
exportButton = uibutton(toolbar,'Text','Export Figure...','ButtonPushedFcn',@exportPlot);
plotGrid = uigridlayout(layout,[1 2],'ColumnWidth',{'1x','1x'},'RowHeight',{'1x'},'Padding',[0 0 0 0]);
uilabel(layout,'Text', ...
    ['Bootstrap: 5,000 whole-eye resamples; angle populations stay separate. Pointwise 95% percentile intervals. ' ...
    'Preview assumes independent eyes; shared wells/batches are not modeled.'], ...
    'WordWrap','on','FontSize',11);
refresh([],[]);
drawnow;
destination = fullfile(root,'comparison_preview','example_groups_ci_gui.png');
exportapp(fig,destination);
writePlot(fullfile(root,'comparison_preview','example_groups_ci_comparison.png'));
fprintf('Saved GUI preview: %s\n',destination);

    function editGroup(~,event)
        included(event.Indices(1)) = logical(event.NewData);
        refresh([],[]);
    end

    function modes = currentModes()
        if strcmp(methodMenu.Value,'both'), modes = {'legacy','bootstrap'};
        else, modes = {methodMenu.Value}; end
    end

    function refresh(~,~)
        delete(plotGrid.Children);
        modes = currentModes(); plotGrid.ColumnWidth = repmat({'1x'},1,numel(modes));
        plotGrid.RowHeight = {'1x'};
        exportButton.Enable = matlab.lang.OnOffSwitchState(any(included));
        for m = 1:numel(modes)
            panel = uipanel(plotGrid,'BorderType','none');
            panel.Layout.Row = 1; panel.Layout.Column = m;
            tile = tiledlayout(panel,1,1,'TileSpacing','compact','Padding','compact');
            ax = nexttile(tile);
            drawProfile(ax,modes{m});
            ax.Tag = ['PreviewAxes_' modes{m}];
        end
    end

    function drawProfile(ax,mode)
        hold(ax,'on'); selected = find(included); curves = gobjects(numel(selected),1);
        for k = 1:numel(selected)
            g = selected(k); group = groups(g);
            if strcmp(mode,'legacy'), bounds = group.Legacy;
            else, bounds = group.Bootstrap; end
            patch(ax,[group.Angle;flipud(group.Angle)], ...
                [bounds(:,1);flipud(bounds(:,2))],colors(g,:), ...
                'FaceAlpha',.16,'EdgeColor','none','HandleVisibility','off', ...
                'Tag',['PreviewBand_' mode],'UserData',group.Name);
        end
        for k = 1:numel(selected)
            g = selected(k); group = groups(g);
            curves(k) = plot(ax,group.Angle,group.Median,'Color',colors(g,:), ...
                'LineWidth',1.9,'DisplayName',sprintf('%s (n = %d)',group.Name,group.Eyes));
        end
        xlim(ax,[0 180]); ylim(ax,[0 255]); grid(ax,'on');
        xlabel(ax,'Angular position (degrees)'); ylabel(ax,'Median intensity (0-255)');
        if strcmp(mode,'legacy')
            title(ax,{'Pixel-based median CI';'Nominal 95%; pixels treated as independent'},'FontSize',12);
        else
            title(ax,{'Whole-eye bootstrap CI';'Pointwise 95%; identical median curves'},'FontSize',12);
        end
        ax.FontSize = 11;
        if ~isempty(curves)
            legend(ax,curves,'Location','southoutside','Interpreter','none','NumColumns',2);
        end
        hold(ax,'off');
    end

    function exportPlot(~,~)
        [name,folder] = uiputfile({'*.png','PNG';'*.pdf','PDF'},'Export CI comparison', ...
            fullfile(root,'comparison_preview','example_groups_ci_comparison.png'));
        if isequal(name,0), return; end
        try
            writePlot(fullfile(folder,name));
        catch problem
            uialert(fig,problem.message,'Export failed');
        end
    end

    function writePlot(destination)
        modes = currentModes();
        output = figure('Visible','off','Color','w','Position',[100 100 1450 670]);
        cleanup = onCleanup(@() close(output));
        tiles = tiledlayout(output,1,numel(modes),'TileSpacing','compact','Padding','compact');
        for m = 1:numel(modes), drawProfile(nexttile(tiles),modes{m}); end
        title(tiles,'Example groups | Pooled-pixel median profiles | 5-degree bins','FontSize',16);
        xlabel(tiles,{'5,000 whole-eye resamples; the same sampled eyes are used across all separate angle populations.', ...
            'Illustrative pointwise intervals: assumes independent eyes; shared wells/batches are not modeled.'},'FontSize',11);
        exportgraphics(output,destination,'Resolution',170);
    end
end

function ci = bootstrapBounds(data,iterations,seed)
raw = data.RAW_data; eyeNames = string(data.IMGname(:)); count = numel(eyeNames);
[known,eyeIndex] = ismember(string(raw.IMG_name),eyeNames);
assert(all(known) && numel(unique(eyeNames)) == count);
assert(count >= 2);
previous = rng; cleanup = onCleanup(@() rng(previous)); rng(seed,'twister');
[~,draws] = bootstrp(iterations,@mean,(1:count)');
drawNumber = repmat(1:iterations,count,1);
weights = accumarray([draws(:),drawNumber(:)],1,[count iterations]);
estimates = nan(36,iterations);
for b = 1:36
    inBin = raw.Angle >= (b-1)*5 & raw.Angle < b*5;
    values = double(raw.Raw_Data(inBin)); ids = eyeIndex(inBin);
    assert(all(values >= 0 & values <= 255 & values == fix(values)));
    histogram = accumarray([values+1,ids],1,[256 count]);
    assert(all(sum(histogram,1) > 0),'Preview requires all eyes to contribute to each bin.');
    assert(sum(histogram(:)) == data.BIN_5_deg_all.Numel(b));
    assert(histogramMedian(sum(histogram,2)) == data.BIN_5_deg_all.Median(b));
    estimates(b,:) = histogramMedian(histogram*weights);
    % Validate compressed-count resampling against actual duplicated pixels.
    for draw = 1:3
        literal = cell(count,1);
        for e = 1:count, literal{e} = values(ids == draws(e,draw)); end
        assert(median(vertcat(literal{:})) == estimates(b,draw));
    end
end
ci = prctile(estimates,[2.5 97.5],2);
assert(all(ci(:,1) <= ci(:,2)) && all(isfinite(ci(:))));
end

function values = histogramMedian(counts)
cumulative = cumsum(counts,1); total = cumulative(end,:);
values = (sum(cumulative < ceil(total/2),1)+sum(cumulative < floor(total/2)+1,1))/2;
end
