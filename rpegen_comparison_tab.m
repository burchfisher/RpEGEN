function api = rpegen_comparison_tab(tab, root, onBusy)
% Group plots and optional permutation tests, independent of extraction.
if nargin < 3, onBusy = []; end
groups = struct('Name',{},'File',{},'Samples',{},'Raw',{},'Bins',{}, ...
    'Include',{},'EyeNames',{},'Bootstrap',{});
busy = false;
loadingDialog = [];
summaries = {};
selected = [];
selectedRows = [];
layout = uigridlayout(tab,[4 1], 'RowHeight',{36,126,36,'1x'});
buttons = uigridlayout(layout,[1 5], 'ColumnWidth',{130,145,145,130,'1x'},'Padding',[0 0 0 0]);
addButton = uibutton(buttons,'Text','Add Groups...', 'ButtonPushedFcn',@addFiles, ...
    'Interruptible','off','BusyAction','cancel','Tag','AddComparisonGroups');
currentButton = uibutton(buttons,'Text','Add Current Results', 'Enable','off', ...
    'Interruptible','off','BusyAction','cancel','Tag','AddCurrentComparisonGroup');
exampleButton = uibutton(buttons,'Text','Example Groups', 'ButtonPushedFcn',@addExamples, 'Tag','ExampleGroups', ...
    'Interruptible','off','BusyAction','cancel');
removeButton = uibutton(buttons,'Text','Remove Selected', 'Enable','off', 'ButtonPushedFcn',@removeGroups);
exportButton = uibutton(buttons,'Text','Export Figure...', 'Enable','off','ButtonPushedFcn',@exportFigure);
groupTable = uitable(layout,'Data',cell(0,4),'ColumnName',{'Show','Group','Animals / eyes','File'}, ...
    'ColumnEditable',[true,true,false,false], 'ColumnFormat',{'logical','char','numeric','char'}, ...
    'ColumnWidth',{55,210,110,'auto'}, 'CellEditCallback',@editGroup, ...
    'CellSelectionCallback',@selectRows,'Tag','ComparisonGroups');
settings = uigridlayout(layout,[1 6], 'ColumnWidth',{70,100,95,100,65,'1x'},'Padding',[0 0 0 0]);
uilabel(settings,'Text','Angle bin');
angleMenu = uidropdown(settings,'Items',{'1 degree','5 degrees','10 degrees'}, ...
    'ItemsData',[1,5,10],'Value',5,'ValueChangedFcn',@refresh,'Tag','HeatmapAngleBin');
uilabel(settings,'Text','Intensity bin');
intensityMenu = uidropdown(settings,'Items',{'1','5','15'},'ItemsData',[1,5,15], ...
    'Value',5,'ValueChangedFcn',@refresh,'Tag','HeatmapIntensityBin');
uilabel(settings,'Text','Color');
normalizationMenu = uidropdown(settings,'Items',{'Pixel count','Percent per angle'}, ...
    'ValueChangedFcn',@refresh,'Tag','HeatmapNormalization');
angleMenu.Tooltip = 'Angular bin for heatmaps and median profiles. Permutation testing has its own bin setting.';
intensityMenu.Tooltip = 'Pixel-intensity interval for heatmaps only; does not change median profiles or permutation tests.';
normalizationMenu.Tooltip = 'Pixel count or percentage within each angle column. Both pool pixels, not equally weighted eyes.';
groupTable.Tooltip = 'Show affects descriptive plots only. Animals / eyes counts saved eye IDs.';
views = uitabgroup(layout,'Tag','ComparisonViews');
heatTab = uitab(views,'Title','Heatmaps');
medianTab = uitab(views,'Title','Median Profiles');
permutationTab = uitab(views,'Title','Permutation Test');
heatGrid = uigridlayout(heatTab,[1 1], 'RowHeight',{'1x'},'Scrollable','on', ...
    'Tag','ComparisonHeatmapGrid');
medianGrid = uigridlayout(medianTab,[3 1],'RowHeight',{30,'1x',36});
medianToolbar = uigridlayout(medianGrid,[1 3],'ColumnWidth',{280,260,'1x'},'Padding',[0 0 0 0]);
confidenceToggle = uicheckbox(medianToolbar,'Text','Pixel-based median CI (nominal 95%)', ...
    'Tooltip','How precisely is the pooled median defined IF pixels are independent? Not the spread of intensities. See ? for an example.', ...
    'Value',true,'ValueChangedFcn',@refreshMedians,'Tag','MedianConfidenceIntervals');
bootstrapToggle = uicheckbox(medianToolbar,'Text','Whole-eye bootstrap CI (95%)', ...
    'Tooltip','How could the same pooled median change with a different sample of eyes? Resamples whole eyes, not separate pixels. See ? for an example.', ...
    'Value',true,'ValueChangedFcn',@refreshMedians,'Tag','BootstrapConfidenceIntervals');
medianPlotGrid = uigridlayout(medianGrid,[1 1],'Padding',[0 0 0 0], ...
    'Tag','MedianPlotGrid');
medianPanel = uipanel(medianPlotGrid,'BorderType','none');
medianTiles = tiledlayout(medianPanel,1,1,'TileSpacing','compact','Padding','compact');
medianAxes = nexttile(medianTiles);
ciStatus = uilabel(medianGrid,'Text', ...
    'Pointwise intervals. Bootstrap assumes independent eyes; shared wells/batches are not modeled.', ...
    'WordWrap','on','Tag','MedianCIStatus');
uilabel(heatGrid,'Text','No groups loaded','HorizontalAlignment','center');
permutation = rpegen_permutation_tab(permutationTab,root,@permutationBusy,@updateExport);
views.SelectionChangedFcn = @changeComparisonView;
api.SetBusy = @setBusy;
api.SetCurrent = @setCurrent;
api.AddFiles = @loadFiles;
api.GetFiles = @getFiles;
api.Export = @writeFigure;
api.Permutation = permutation;
renderMedians(false);

    function files = getFiles()
        files = {groups.File};
    end

    function setCurrent(file)
        currentButton.ButtonPushedFcn = @(~,~) addCurrent(file);
        currentButton.Enable = matlab.lang.OnOffSwitchState(~busy);
    end

    function setBusy(value)
        busy = value;
        permutation.SetBusy(value);
        controls = {addButton,exampleButton,angleMenu,intensityMenu,normalizationMenu,confidenceToggle,bootstrapToggle};
        for k = 1:numel(controls)
            controls{k}.Enable = matlab.lang.OnOffSwitchState(~value);
        end
        currentButton.Enable = matlab.lang.OnOffSwitchState(~value && ~isempty(currentButton.ButtonPushedFcn));
        removeButton.Enable = matlab.lang.OnOffSwitchState(~value && ~isempty(selectedRows) && views.SelectedTab~=permutationTab);
        updateExport();
        if value
            groupTable.Enable = 'off';
        else
            groupTable.Enable = 'on';
        end
    end

    function permutationBusy(value)
        if isempty(onBusy), setBusy(value); else, onBusy(value); end
    end

    function updateExport()
        ready = ~isempty(selected);
        if views.SelectedTab==permutationTab, ready = permutation.HasResults(); end
        exportButton.Enable = matlab.lang.OnOffSwitchState(~busy && ready);
    end

    function changeComparisonView(~,~)
        if views.SelectedTab==permutationTab
            settings.Visible = 'off'; groupTable.Visible = 'off'; layout.RowHeight = {36,0,0,'1x'};
        else
            settings.Visible = 'on'; groupTable.Visible = 'on'; layout.RowHeight = {36,126,36,'1x'};
        end
        removeButton.Enable = matlab.lang.OnOffSwitchState(~busy && ~isempty(selectedRows) && views.SelectedTab~=permutationTab);
        updateExport();
    end

    function addFiles(~,~)
        if busy, return; end
        [names,folder] = uigetfile('*.mat','Select group results',root,'MultiSelect','on');
        if isequal(names,0), return; end
        if ischar(names), names = {names}; end
        try
            loadFiles(fullfile(folder,names),addButton);
        catch problem
            alert(problem);
        end
    end

    function addCurrent(file)
        if busy, return; end
        try
            loadFiles({file},currentButton);
        catch problem
            alert(problem);
        end
    end

    function addExamples(~,~)
        if busy, return; end
        names = {'DMSO_9dpf.mat','IWR1_9dpf.mat','DMSO_4dpi.mat','IWR1_4dpi.mat'};
        try
            loadFiles(fullfile(root,'JoVE_dataset','Output',names),exampleButton);
        catch problem
            alert(problem);
        end
    end

    function loadFiles(files,loadingButton)
        if busy, return; end
        if nargin<2, loadingButton = []; end
        files = cellstr(string(files));
        for k = 1:numel(files)
            files{k} = char(java.io.File(files{k}).getCanonicalPath());
        end
        files = unique(files,'stable');
        files = files(~ismember(files,{groups.File}));
        if isempty(files), return; end
        buttonText = '';
        if ~isempty(loadingButton), buttonText = loadingButton.Text; loadingButton.Text = 'Loading...'; end
        setBusy(true);
        cleanup = onCleanup(@() finishLoading(loadingButton,buttonText));
        parentFigure = ancestor(tab,'figure');
        if strcmp(parentFigure.Visible,'on')
            loadingDialog = uiprogressdlg(parentFigure,'Title','Loading Groups', ...
                'Message','Preparing group files...','Indeterminate','on','Cancelable','off');
        end
        drawnow nocallbacks;
        if ~isempty(onBusy), onBusy(true); end
        pending = groups;
        for k = 1:numel(files)
            file = files{k};
            [~,name,extension] = fileparts(file);
            loadingMessage(sprintf('Loading group %d of %d: %s%s',k,numel(files),name,extension));
            loaded = load(file);
            names = fieldnames(loaded); candidate = []; groupName = '';
            required = {'RAW_data','IMGname','BIN_1_deg_all','BIN_5_deg_all','BIN_10_deg_all'};
            for j = 1:numel(names)
                entry = loaded.(names{j});
                if isstruct(entry) && isscalar(entry) && all(isfield(entry,required))
                    candidate = entry; groupName = names{j}; break;
                end
            end
            if isempty(candidate)
                error('RpEGEN:Comparison','%s is not a group result MAT file.',files{k});
            end
            rpegen_heatmap(candidate.RAW_data,5,5,'Pixel count');
            bins = {candidate.BIN_1_deg_all,candidate.BIN_5_deg_all,candidate.BIN_10_deg_all};
            for j = 1:3
                if ~istable(bins{j}) || ~all(ismember({'Angle','Median','CI95L','CI95U'},bins{j}.Properties.VariableNames))
                    error('RpEGEN:Comparison','%s must contain Angle, Median, CI95L, and CI95U in each binned table.',files{k});
                end
            end
            pending(end+1) = struct('Name',groupName,'File',file, ...
                'Samples',numel(candidate.IMGname),'Raw',candidate.RAW_data, ...
                'Bins',{bins},'Include',true,'EyeNames',{candidate.IMGname}, ...
                'Bootstrap',{{[],[],[]}}); %#ok<AGROW>
        end
        loadingMessage('Preparing heatmaps and median profiles...');
        groups = pending; updateTable(); refresh([],[]);
    end

    function loadingMessage(message)
        if ~isempty(loadingDialog) && isvalid(loadingDialog), loadingDialog.Message = message; end
        drawnow nocallbacks;
    end

    function finishLoading(loadingButton,buttonText)
        if ~isempty(loadingDialog) && isvalid(loadingDialog), close(loadingDialog); end
        loadingDialog = [];
        if ~isempty(loadingButton) && isvalid(loadingButton), loadingButton.Text = buttonText; end
        setBusy(false);
        if ~isempty(onBusy), onBusy(false); end
    end

    function updateTable()
        rows = cell(numel(groups),4);
        for k = 1:numel(groups)
            rows(k,:) = {groups(k).Include,groups(k).Name,groups(k).Samples,groups(k).File};
        end
        groupTable.Data = rows;
        selectedRows = []; removeButton.Enable = 'off';
    end

    function selectRows(~,event)
        selectedRows = unique(event.Indices(:,1));
        removeButton.Enable = matlab.lang.OnOffSwitchState(~busy && ~isempty(selectedRows));
    end

    function editGroup(~,event)
        if busy, return; end
        row = event.Indices(1); column = event.Indices(2);
        if column == 1
            groups(row).Include = logical(event.NewData);
        else
            name = strtrim(char(string(event.NewData)));
            if isempty(name)
                groupTable.Data{row,2} = groups(row).Name;
                return;
            end
            groups(row).Name = name;
            groupTable.Data{row,2} = name;
        end
        refresh([],[]);
    end

    function removeGroups(~,~)
        if busy || isempty(selectedRows), return; end
        groups(selectedRows) = [];
        updateTable(); refresh([],[]);
    end

    function refresh(~,~)
        selected = find([groups.Include]);
        permutation.SetGroups(groups);
        delete(heatGrid.Children);
        summaries = cell(1,numel(selected));
        updateExport();
        if isempty(selected)
            heatGrid.RowHeight = {'1x'}; heatGrid.ColumnWidth = {'1x'};
            uilabel(heatGrid,'Text','No groups selected','HorizontalAlignment','center');
            renderMedians(false);
            return;
        end
        columns = min(2,numel(selected)); rows = ceil(numel(selected)/columns);
        heatGrid.ColumnWidth = repmat({'1x'},1,columns);
        if rows == 1, heatGrid.RowHeight = {'1x'};
        else, heatGrid.RowHeight = repmat({310},1,rows); end
        for k = 1:numel(selected)
            summaries{k} = rpegen_heatmap(groups(selected(k)).Raw,angleMenu.Value, ...
                intensityMenu.Value,normalizationMenu.Value);
        end
        limit = sharedLimit();
        for k = 1:numel(selected)
            panel = uipanel(heatGrid,'BorderType','none');
            panel.Layout.Row = ceil(k/columns); panel.Layout.Column = mod(k-1,columns)+1;
            tile = tiledlayout(panel,1,1,'TileSpacing','compact','Padding','compact');
            ax = nexttile(tile);
            drawHeatmap(ax,k,limit);
            ax.Tag = 'GroupHeatmap';
        end
        scroll(heatGrid,'top');
        renderMedians(false);
    end

    function limit = sharedLimit()
        limit = max(1,max(cellfun(@(s) max(s.Values(:)),summaries)));
    end

    function bar = drawHeatmap(ax,k,limit)
        summary = summaries{k}; group = groups(selected(k));
        x = (summary.AngleEdges(1:end-1)+summary.AngleEdges(2:end))/2;
        y = (summary.IntensityEdges(1:end-1)+summary.IntensityEdges(2:end))/2;
        imagesc(ax,x,y,summary.Values); ax.YDir = 'normal';
        xlim(ax,[0 180]); ylim(ax,[0 255]); clim(ax,[0 limit]); colormap(ax,turbo(256));
        bar = colorbar(ax); bar.Label.String = summary.ColorLabel;
        title(ax,sprintf('%s (n = %d)',group.Name,group.Samples),'Interpreter','none');
        xlabel(ax,'Angular distance (degrees)'); ylabel(ax,'Pixel intensity');
    end

    function drawMedians(ax)
        % cla alone preserves bands whose handles are hidden from the legend.
        delete(findall(ax,'Tag','MedianConfidenceBand'));
        delete(findall(ax,'Tag','BootstrapConfidenceBand'));
        cla(ax); legend(ax,'off'); hold(ax,'on'); colors = lines(numel(groups));
        index = find([1,5,10] == angleMenu.Value,1);
        % Broad bootstrap bands first, darker legacy bands next, medians last.
        for method = {'bootstrap','legacy'}
            mode = method{1};
            if strcmp(mode,'bootstrap')
                if ~bootstrapToggle.Value, continue; end
                opacity = 0.16; bandTag = 'BootstrapConfidenceBand';
            else
                if ~confidenceToggle.Value, continue; end
                opacity = 0.45; bandTag = 'MedianConfidenceBand';
            end
            for k = 1:numel(selected)
                bins = groups(selected(k)).Bins{index};
                x = bins.Angle(:); lower = bins.CI95L(:); upper = bins.CI95U(:);
                if strcmp(mode,'bootstrap')
                    ci = groups(selected(k)).Bootstrap{index};
                    if isempty(ci) || isfield(ci,'Error'), continue; end
                    lower = ci.Lower; upper = ci.Upper;
                end
                valid = isfinite(x) & isfinite(lower) & isfinite(upper) & lower <= upper;
                edges = diff([false; valid; false]);
                starts = find(edges == 1); stops = find(edges == -1)-1;
                for segment = 1:numel(starts)
                    rows = starts(segment):stops(segment);
                    patch(ax,[x(rows);flipud(x(rows))], ...
                        [lower(rows);flipud(upper(rows))],colors(selected(k),:), ...
                        'FaceAlpha',opacity,'EdgeColor','none','HandleVisibility','off', ...
                        'Tag',bandTag,'UserData',groups(selected(k)).Name);
                end
            end
        end
        for k = 1:numel(selected)
            group = groups(selected(k)); bins = group.Bins{index};
            plot(ax,bins.Angle,bins.Median,'LineWidth',1.6,'Color',colors(selected(k),:), ...
                'DisplayName',sprintf('%s (n = %d)',group.Name,group.Samples));
        end
        hold(ax,'off'); xlim(ax,[0 180]); ylim(ax,[0 255]); grid(ax,'on');
        title(ax,sprintf('Pooled-pixel median profiles (%d-degree bins)',angleMenu.Value),'FontSize',12);
        xlabel(ax,'Angular distance (degrees)'); ylabel(ax,'Median intensity');
        if ~isempty(selected)
            legend(ax,'Location','southoutside','Interpreter','none','NumColumns',2);
        end
    end

    function refreshMedians(~,~)
        if busy, return; end
        renderMedians(true);
    end

    function prepareBootstrap()
        if ~bootstrapToggle.Value, return; end
        index = find([1 5 10]==angleMenu.Value,1);
        missing = selected(arrayfun(@(k) isempty(groups(k).Bootstrap{index}),selected));
        if isempty(missing), return; end
        wasBusy = busy; setBusy(true);
        cleanup = onCleanup(@() setBusy(wasBusy));
        for k = missing
            ciStatus.Text = sprintf('Bootstrapping %s: 5,000 whole-eye resamples...',groups(k).Name);
            if ~isempty(loadingDialog) && isvalid(loadingDialog)
                loadingDialog.Message = sprintf('Calculating confidence bands: %s (5,000 whole-eye resamples)',groups(k).Name);
            end
            drawnow nocallbacks;
            try
                ci = rpegen_bootstrap_ci(groups(k).Raw,groups(k).EyeNames,angleMenu.Value);
                bins = groups(k).Bins{index};
                if ~isequaln(ci.Angle,bins.Angle(:)) || ~isequaln(ci.Median,double(bins.Median(:)))
                    error('RpEGEN:Bootstrap','Saved bins do not match the raw pooled medians.');
                end
                groups(k).Bootstrap{index} = ci;
            catch problem
                groups(k).Bootstrap{index} = struct('Error',problem.message);
            end
        end
    end

    function renderMedians(preserveView)
        limits = [medianAxes.XLim;medianAxes.YLim];
        prepareBootstrap();
        drawMedians(medianAxes); medianAxes.Tag = 'GroupMedianAxes';
        if preserveView
            xlim(medianAxes,limits(1,:)); ylim(medianAxes,limits(2,:));
        end
        drawnow nocallbacks;
        if bootstrapToggle.Value
            ciStatus.Text = ['Pointwise intervals; 5,000 resamples. Bootstrap assumes independent eyes; ' ...
                'shared wells/batches are not modeled.'];
        elseif confidenceToggle.Value
            ciStatus.Text = 'Pixel-based median intervals: nominal 95%, pointwise; not between-animal confidence intervals.';
        else
            ciStatus.Text = 'Pooled-pixel medians; confidence intervals hidden.';
        end
        issues = {};
        if bootstrapToggle.Value
            index = find([1 5 10]==angleMenu.Value,1);
            for k = selected
                ci = groups(k).Bootstrap{index};
                if isfield(ci,'Error')
                    issues{end+1} = sprintf('%s: %s',groups(k).Name,ci.Error); %#ok<AGROW>
                elseif any(~isfinite(ci.Lower))
                    issues{end+1} = sprintf('%s: %d bins without a bootstrap CI', ...
                        groups(k).Name,nnz(~isfinite(ci.Lower))); %#ok<AGROW>
                end
            end
        end
        if ~isempty(issues), ciStatus.Text = strjoin(issues,'; '); end
        ciStatus.Tooltip = ciStatus.Text;
    end

    function exportFigure(~,~)
        if busy, return; end
        if views.SelectedTab==permutationTab
            if ~permutation.HasResults(), return; end
            destination = fileparts(groups(1).File);
        else
            if isempty(selected), return; end
            destination = fileparts(groups(selected(1)).File);
        end
        [name,folder] = uiputfile({'*.pdf','Vector PDF'},'Export comparison', ...
            fullfile(destination, ...
            [lower(strrep(views.SelectedTab.Title,' ','_')) '.pdf']));
        if isequal(name,0), return; end
        try
            writeFigure(fullfile(folder,name),views.SelectedTab.Title);
        catch problem
            alert(problem);
        end
    end

    function writeFigure(file,view)
        if strcmp(view,'Permutation Test')
            permutation.ExportFigure(file); return;
        end
        if isempty(selected), error('RpEGEN:Comparison','Select at least one group.'); end
        count = numel(selected); rows = ceil(count/2);
        exportFig = figure('Visible','off','Color','w','Position',[100 100 1100 max(450,330*rows)]);
        cleanup = onCleanup(@() close(exportFig));
        if strcmp(view,'Heatmaps')
            tiles = tiledlayout(exportFig,rows,min(2,count),'TileSpacing','compact');
            for k = 1:count, drawHeatmap(nexttile(tiles),k,sharedLimit()); end
        elseif strcmp(view,'Median Profiles')
            prepareBootstrap();
            exportFig.Position(4) = 600;
            tiles = tiledlayout(exportFig,1,1,'TileSpacing','compact');
            drawMedians(nexttile(tiles));
            methods = {};
            if confidenceToggle.Value, methods{end+1} = 'Pixel-based median CI (nominal 95%, dark)'; end
            if bootstrapToggle.Value, methods{end+1} = 'Whole-eye bootstrap CI (95%, light)'; end
            if isempty(methods), methods = {'Confidence intervals hidden'}; end
            xlabel(tiles,{strjoin(methods,'; ');ciStatus.Text},'FontSize',9,'Interpreter','none');
        else
            error('RpEGEN:Comparison','Unknown comparison view.');
        end
        rpegen_export_figure(exportFig,file);
    end

    function alert(problem)
        uialert(ancestor(tab,'figure'),problem.message,'Group Comparison');
    end
end
