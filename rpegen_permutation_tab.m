function api = rpegen_permutation_tab(tab,root,onBusy,onStateChanged)
% Interactive adapter for the unchanged pixel-level permutation test.
groups = struct('Name',{},'File',{},'Samples',{},'Raw',{});
running = false; blocked = false; cancelRequested = false;
report = table; metadata = struct; started = [];
layout = uigridlayout(tab,[7 1],'RowHeight',{32,32,32,28,24,'1x',36});
pair = uigridlayout(layout,[1 4],'ColumnWidth',{60,'1x',60,'1x'},'Padding',[0 0 0 0]);
uilabel(pair,'Text','Group A');
groupA = uidropdown(pair,'Items',{'No groups loaded'},'ItemsData',{''}, ...
    'ValueChangedFcn',@settingsChanged,'Tag','PermutationGroupA');
uilabel(pair,'Text','Group B');
groupB = uidropdown(pair,'Items',{'No groups loaded'},'ItemsData',{''}, ...
    'ValueChangedFcn',@settingsChanged,'Tag','PermutationGroupB');
settings = uigridlayout(layout,[1 8],'ColumnWidth',{62,95,44,95,90,85,100,'1x'},'Padding',[0 0 0 0]);
uilabel(settings,'Text','Angle bin');
angle = uidropdown(settings,'Items',{'1 degree','5 degrees','10 degrees'}, ...
    'ItemsData',[1 5 10],'Value',1,'ValueChangedFcn',@settingsChanged,'Tag','PermutationBin');
uilabel(settings,'Text','Mode');
mode = uidropdown(settings,'Items',{'Auto','Manual'},'Value','Auto', ...
    'ValueChangedFcn',@settingsChanged,'Tag','PermutationMode');
uilabel(settings,'Text','Budget (min)');
budget = uispinner(settings,'Value',5,'Limits',[0.1 60],'Step',1, ...
    'Tooltip','Auto aims for 45,618 shuffles per bin; slower bins may use fewer to fit this approximate total runtime, with a minimum of 10,000. See ?.', ...
    'ValueChangedFcn',@settingsChanged,'Tag','PermutationBudget');
uilabel(settings,'Text','Manual draws');
iterations = uispinner(settings,'Value',5000,'Limits',[100 1e7],'Step',1000, ...
    'Tooltip','One draw is one random pixel-label shuffle within an angular bin. Manual uses this exact count per bin; it is not an animal count.', ...
    'RoundFractionalValues','on','ValueChangedFcn',@settingsChanged,'Tag','PermutationIterations');
actions = uigridlayout(layout,[1 5],'ColumnWidth',{105,85,160,160,'1x'},'Padding',[0 0 0 0]);
runButton = uibutton(actions,'Text','Run Test','ButtonPushedFcn',@runClicked,'Tag','RunPermutation');
cancelButton = uibutton(actions,'Text','Cancel','ButtonPushedFcn',@cancelRun,'Tag','CancelPermutation');
parallelToggle = uicheckbox(actions,'Text','Parallel CPU (up to 4)', ...
    'Value',logical(license('test','Distrib_Computing_Toolbox')), ...
    'ValueChangedFcn',@settingsChanged,'Tag','PermutationParallel');
exportButton = uibutton(actions,'Text','Export Results...', ...
    'ButtonPushedFcn',@exportResults,'Tag','ExportPermutationResults');
displayOptions = uigridlayout(layout,[1 4],'ColumnWidth',{170,190,170,'1x'},'Padding',[0 0 0 0]);
focusToggle = uicheckbox(displayOptions,'Text','Focus near p = 0.05','Value',false, ...
    'ValueChangedFcn',@changeFocus,'Tag','PermutationFocus');
intervalToggle = uicheckbox(displayOptions,'Text','Show simulation intervals','Value',false, ...
    'Tooltip','Uncertainty in the p-value from a finite number of random shuffles, not uncertainty from sampling different eyes. See ? for an example.', ...
    'ValueChangedFcn',@changeDisplay,'Tag','PermutationIntervals');
shadeToggle = uicheckbox(displayOptions,'Text','Shade p < 0.05','Value',true, ...
    'Tooltip','Highlights unadjusted pixel-level p-values below 0.05. Does not show effect size or correct for testing many bins.', ...
    'ValueChangedFcn',@changeDisplay,'Tag','PermutationShading');
status = uilabel(layout,'Text','Load at least two groups.','Tag','PermutationStatus');
views = uitabgroup(layout,'Tag','PermutationViews');
plotTab = uitab(views,'Title','P-value Plot'); dataTab = uitab(views,'Title','Results Table');
plotGrid = uigridlayout(plotTab,[1 1],'Padding',[0 0 0 0]);
panel = uipanel(plotGrid,'BorderType','none');
tiles = tiledlayout(panel,1,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(tiles);
dataGrid = uigridlayout(dataTab,[1 1],'Padding',[0 0 0 0]);
resultsTable = uitable(dataGrid,'Data',table,'ColumnEditable',false,'Tag','PermutationResults');
resultsTable.Tooltip = 'LegacyP = Exceedances / Iterations. MC95L/U describe simulation precision only; NearThreshold means the interval includes 0.05. See ? for definitions.';
note = ['Pixel-level test; unadjusted p-values. Optional bars: pointwise 95% Monte Carlo intervals, ' ...
    'not biological CIs. Zero = no exceedances observed.'];
uilabel(layout,'Text',note,'WordWrap','on');
api.SetGroups = @setGroups; api.SetBusy = @setBusy; api.Run = @runTest;
api.Cancel = @cancelRun; api.HasResults = @hasResults;
api.ExportResults = @writeResults; api.ExportFigure = @writeFigure;
api.GetResults = @getResults;
drawPlot(ax); updateControls();

    function setGroups(incoming)
        if running, return; end
        names = {incoming.Name}; files = {incoming.File};
        if isequal(names,{groups.Name}) && isequal(files,{groups.File}), return; end
        previousA = groupA.Value; previousB = groupB.Value; groups = incoming;
        if isempty(groups)
            labels = {'No groups loaded'}; files = {''};
        else
            labels = arrayfun(@(k) sprintf('%d. %s (n = %d)',k,groups(k).Name, ...
                groups(k).Samples),1:numel(groups),'UniformOutput',false);
        end
        groupA.ItemsData = {}; groupB.ItemsData = {};
        groupA.Items = labels; groupB.Items = labels;
        groupA.ItemsData = files; groupB.ItemsData = files;
        if ismember(previousA,files), groupA.Value = previousA; else, groupA.Value = files{1}; end
        if ismember(previousB,files), groupB.Value = previousB;
        else, groupB.Value = files{min(2,numel(files))}; end
        if numel(files)>1 && strcmp(groupA.Value,groupB.Value)
            groupB.Value = files{find(~strcmp(files,groupA.Value),1)};
        end
        settingsChanged([],[]);
    end

    function setBusy(value)
        blocked = value; updateControls();
    end

    function updateControls()
        available = ~running && ~blocked;
        controls = {groupA,groupB,angle,mode,parallelToggle};
        for k = 1:numel(controls), controls{k}.Enable = matlab.lang.OnOffSwitchState(available); end
        budget.Enable = matlab.lang.OnOffSwitchState(available && strcmp(mode.Value,'Auto'));
        iterations.Enable = matlab.lang.OnOffSwitchState(available && strcmp(mode.Value,'Manual'));
        runButton.Enable = matlab.lang.OnOffSwitchState(available && numel(groups)>=2 && ~strcmp(groupA.Value,groupB.Value));
        cancelButton.Enable = matlab.lang.OnOffSwitchState(running && ~cancelRequested);
        focusToggle.Enable = matlab.lang.OnOffSwitchState(available);
        intervalToggle.Enable = matlab.lang.OnOffSwitchState(available);
        shadeToggle.Enable = matlab.lang.OnOffSwitchState(available);
        exportButton.Enable = matlab.lang.OnOffSwitchState(available && ~isempty(report));
    end

    function settingsChanged(~,~)
        if running, return; end
        report = table; metadata = struct; resultsTable.Data = table;
        drawPlot(ax); updateControls(); onStateChanged();
        if numel(groups)<2
            status.Text = 'Load at least two groups.';
        elseif strcmp(groupA.Value,groupB.Value)
            status.Text = 'Select two different groups.';
        else
            status.Text = 'Ready';
        end
    end

    function runClicked(~,~)
        try
            runTest();
        catch problem
            if ~strcmp(problem.identifier,'RpEGEN:Cancelled') && isvalid(tab)
                uialert(ancestor(tab,'figure'),problem.message,'Permutation Test');
            end
        end
    end

    function runTest()
        if running || blocked, return; end
        assert(numel(groups)>=2 && ~strcmp(groupA.Value,groupB.Value), ...
            'RpEGEN:Permutation','Select two different groups.');
        a = groups(strcmp({groups.File},groupA.Value)); b = groups(strcmp({groups.File},groupB.Value));
        settingsChanged([],[]); cancelRequested = false; running = true; started = tic;
        cleanup = onCleanup(@() finishRun());
        onBusy(true); updateControls();
        if strcmp(mode.Value,'Auto'), reps = 'auto'; else, reps = iterations.Value; end
        options = struct('TimeBudgetSeconds',60*budget.Value,'UseParallel',parallelToggle.Value, ...
            'NumWorkers',4,'ProgressFcn',@progress,'CancelFcn',@cancelled);
        try
            [newReport,runInfo] = rpegen_permutation_profile(a.Raw,b.Raw,angle.Value,reps,options);
            if cancelled(), error('RpEGEN:Cancelled','Permutation test cancelled.'); end
            report = newReport;
            metadata = struct('GroupA',a.Name,'GroupB',b.Name,'FileA',a.File,'FileB',b.File, ...
                'EyesA',a.Samples,'EyesB',b.Samples,'Mode',mode.Value,'RequestedIterations',reps, ...
                'RunInfo',runInfo,'CompletedUTC',char(datetime('now','TimeZone','UTC')), ...
                'AutoMinIterations',10000,'AutoMaxIterations',50000,'AutoHalfWidth',0.002, ...
                'PermutationUnit','pixel','PValueFormula','Exceedances / Iterations', ...
                'Threshold',0.05,'MultipleTestingAdjustment','none');
            resultsTable.Data = report; drawPlot(ax);
            status.Text = sprintf('Complete: %d bins | %.1f s | %d CPU workers | %d-%d draws/bin | %d uncertain at 0.05 | %d zeros', ...
                height(report),toc(started),runInfo.Workers,min(report.Iterations),max(report.Iterations), ...
                nnz(report.NearThreshold),nnz(report.Exceedances==0));
            status.Tooltip = status.Text;
        catch problem
            if isvalid(tab)
                if strcmp(problem.identifier,'RpEGEN:Cancelled'), status.Text = 'Cancelled. No partial results retained.';
                else, status.Text = ['Failed: ' problem.message]; end
                status.Tooltip = status.Text;
            end
            rethrow(problem);
        end
    end

    function progress(done,total,message)
        if ~isvalid(tab), return; end
        status.Text = sprintf('%s | %d / %d bins (%d%%) | %.0f s',message,done,total,round(100*done/total),toc(started));
        drawnow limitrate;
    end

    function value = cancelled()
        drawnow limitrate;
        value = cancelRequested || ~isvalid(tab);
    end

    function cancelRun(varargin)
        if ~running, return; end
        cancelRequested = true; status.Text = 'Cancelling...'; updateControls();
    end

    function finishRun()
        running = false;
        if isvalid(tab), onBusy(false); updateControls(); onStateChanged(); end
    end

    function changeFocus(~,~)
        drawPlot(ax);
    end

    function changeDisplay(~,~)
        limitsX = ax.XLim; limitsY = ax.YLim;
        drawPlot(ax);
        xlim(ax,limitsX); ylim(ax,limitsY);
    end

    function drawPlot(target)
        cla(target); hold(target,'on'); legend(target,'off');
        if ~isempty(report)
            edges = [report.Angle-metadata.RunInfo.BinWidth/2; ...
                report.Angle(end)+metadata.RunInfo.BinWidth/2];
            if shadeToggle.Value && any(report.LegacyP<0.05)
                transitions = diff([false; report.LegacyP<0.05; false]);
                left = edges(transitions==1);
                right = edges(transitions==-1);
                patch(target,'XData',[left left right right]', ...
                    'YData',repmat([0;1;1;0],1,numel(left)), ...
                    'FaceColor',[0.1 0.6 0.5],'FaceAlpha',0.12,'EdgeColor','none', ...
                    'HitTest','off','PickableParts','none', ...
                    'DisplayName','p < 0.05 (unadjusted)','Tag','PermutationRibbon');
            end
            stairs(target,edges,[report.LegacyP; report.LegacyP(end)], ...
                'LineWidth',1.25,'Color',[0.05 0.4 0.6], ...
                'DisplayName','Permutation p-value','Tag','PermutationP');
            if intervalToggle.Value
                errorbar(target,report.Angle,report.LegacyP,report.LegacyP-report.MC95L, ...
                    report.MC95U-report.LegacyP,'LineStyle','none','Color',[0.25 0.25 0.25], ...
                    'LineWidth',1.25,'CapSize',5, ...
                    'DisplayName','95% simulation interval','Tag','PermutationMC');
            end
            near = report.NearThreshold;
            if any(near)
                scatter(target,report.Angle(near),report.LegacyP(near),35,[0.85 0.3 0.1], ...
                    'LineWidth',1.25,'MarkerFaceColor','none', ...
                    'DisplayName','Uncertain at 0.05','Tag','PermutationUncertain');
            end
            title(target,{sprintf('%s vs %s',metadata.GroupA,metadata.GroupB); ...
                sprintf('%d-degree bins | n = %d vs %d eyes',metadata.RunInfo.BinWidth,metadata.EyesA,metadata.EyesB)}, ...
                'Interpreter','none','FontSize',12);
            legend(target,'Location','southoutside','NumColumns',3);
        else
            title(target,'Pixel permutation test');
        end
        yline(target,0.05,'--','0.05','HandleVisibility','off');
        xlim(target,[0 180]);
        if focusToggle.Value, ylim(target,[0 0.12]); else, ylim(target,[0 1]); end
        xlabel(target,'Angular distance (degrees)'); ylabel(target,'p-value'); grid(target,'on');
        hold(target,'off'); target.Tag = 'PermutationAxes';
    end

    function exportResults(~,~)
        [name,folder] = uiputfile({'*.csv','CSV';'*.mat','MAT file'},'Export permutation results', ...
            fullfile(root,'permutation_results.csv'));
        if isequal(name,0), return; end
        try
            writeResults(fullfile(folder,name));
        catch problem
            uialert(ancestor(tab,'figure'),problem.message,'Export Results');
        end
    end

    function writeResults(file)
        assert(~isempty(report) && ~running,'RpEGEN:Permutation','No completed test to export.');
        [~,~,extension] = fileparts(file);
        if strcmpi(extension,'.mat')
            save(file,'report','metadata');
        elseif strcmpi(extension,'.csv')
            output = report; count = height(report);
            output.GroupA = repmat(string(metadata.GroupA),count,1);
            output.GroupB = repmat(string(metadata.GroupB),count,1);
            output.FileA = repmat(string(metadata.FileA),count,1);
            output.FileB = repmat(string(metadata.FileB),count,1);
            output.BinWidth = repmat(metadata.RunInfo.BinWidth,count,1);
            output.Seed = repmat(metadata.RunInfo.Seed,count,1);
            output.Mode = repmat(string(metadata.Mode),count,1);
            output.EyesA = repmat(metadata.EyesA,count,1);
            output.EyesB = repmat(metadata.EyesB,count,1);
            output.Workers = repmat(metadata.RunInfo.Workers,count,1);
            output.Execution = repmat(string(metadata.RunInfo.Execution),count,1);
            output.CompletedUTC = repmat(string(metadata.CompletedUTC),count,1);
            output.TimeBudgetSeconds = repmat(metadata.RunInfo.TimeBudgetSeconds,count,1);
            output.PermutationUnit = repmat("pixel",count,1);
            output.MultipleTestingAdjustment = repmat("none",count,1);
            writetable(output,file);
        else
            error('RpEGEN:Permutation','Choose CSV or MAT format.');
        end
    end

    function writeFigure(file)
        assert(~isempty(report) && ~running,'RpEGEN:Permutation','No completed test to export.');
        fig = figure('Visible','off','Color','w','Position',[100 100 1100 600]);
        cleanup = onCleanup(@() close(fig));
        grid = tiledlayout(fig,1,1,'TileSpacing','compact'); drawPlot(nexttile(grid));
        xlabel(grid,note,'FontSize',9); rpegen_export_figure(fig,file);
    end

    function [value,details] = getResults()
        value = report; details = metadata;
    end

    function value = hasResults()
        value = ~isempty(report);
    end
end
