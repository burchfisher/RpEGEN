function fig = RpEGEN_GUI(varargin)
% Launch the setup and quality-control interface for the published pipeline.
root = fileparts(mfilename('fullpath'));
options = inputParser;
addParameter(options, 'Visible', 'on');
addParameter(options, 'ResultsFile', '');
parse(options, varargin{:});
addpath(root, fullfile(root, 'ReadImageJROI-master'));
data = [];
resultPath = '';
qcFolder = '';
busy = false;
helpWindow = [];
fig = uifigure('Name', 'RpEGEN', 'Position', [100 100 1180 780], ...
    'CloseRequestFcn', @closeApp, 'WindowKeyPressFcn',@sampleKey, ...
    'Tag', 'RpEGENGUI', 'Visible', options.Results.Visible,'DeleteFcn',@(~,~) closeHelp());
layout = uigridlayout(fig, [2 1], 'RowHeight', {48, '1x'}, 'Padding', [16 12 16 12]);
header = uigridlayout(layout, [1 3], 'ColumnWidth', {'1x', 300,40}, 'Padding', [0 0 0 0]);
uilabel(header, 'Text', 'RpEGEN', 'FontSize', 24, 'FontWeight', 'bold');
environment = uilabel(header, 'Text', ['MATLAB ' version('-release')], ...
    'HorizontalAlignment', 'right');
uibutton(header,'Text','?','FontSize',20,'Tag','ContextHelp', ...
    'Tooltip','Help for the current tab: controls, plots, and table columns','ButtonPushedFcn',@showHelp);
tabs = uitabgroup(layout);
setupTab = uitab(tabs, 'Title', 'Setup');
qcTab = uitab(tabs, 'Title', 'Quality Control');
resultsTab = uitab(tabs, 'Title', 'Results');
compareTab = uitab(tabs,'Title','Compare Groups');
comparison = rpegen_comparison_tab(compareTab,root,@comparisonBusy);
setup = uigridlayout(setupTab, [8 3], 'ColumnWidth', {150, '1x', 110}, ...
    'RowHeight', {34,34,34,34,38,40,'1x',30});
roiField = folderRow(1, 'ROI folder');
imageField = folderRow(2, 'Image folder');
outputField = folderRow(3, 'Output folder');
roiField.Tag = 'ROIdir'; imageField.Tag = 'IMGdir'; outputField.Tag = 'Outputdir';
label = uilabel(setup, 'Text', 'Group name'); label.Layout.Row = 4; label.Layout.Column = 1;
groupField = uieditfield(setup, 'text', 'Value', 'DMSO_4dpi', 'Tag', 'GroupName');
groupField.Layout.Row = 4; groupField.Layout.Column = 2;
demoButton = uibutton(setup, 'Text', 'Example Data', 'ButtonPushedFcn', @useDemo);
demoButton.Layout.Row = 4; demoButton.Layout.Column = 3;
label = uilabel(setup, 'Text', 'Brightfield page'); label.Layout.Row = 5; label.Layout.Column = 1;
channel = uispinner(setup, 'Limits', [1 Inf], 'Step', 1, 'Value', 3, ...
    'RoundFractionalValues', 'on', 'Tag', 'Channel');
channel.Layout.Row = 5; channel.Layout.Column = 2;
actions = uigridlayout(setup, [1 3], 'ColumnWidth', {150,150,'1x'}, 'Padding', [0 0 0 0]);
actions.Layout.Row = 6; actions.Layout.Column = [1 3];
checkButton = uibutton(actions, 'Text', 'Check Inputs', 'ButtonPushedFcn', @checkInputs);
runButton = uibutton(actions, 'Text', 'Run Analysis', 'ButtonPushedFcn', @runAnalysis);
loadButton = uibutton(actions, 'Text', 'Open Existing Results', 'ButtonPushedFcn', @loadResults);
sampleTable = uitable(setup, 'ColumnName', {'Sample', 'Height', 'Width'}, ...
    'ColumnEditable', false, 'Tag', 'Samples');
sampleTable.Layout.Row = 7; sampleTable.Layout.Column = [1 3];
status = uilabel(setup, 'Text', 'Ready', 'Tag', 'Status');
status.Layout.Row = 8; status.Layout.Column = [1 3];

qc = uigridlayout(qcTab, [3 3], 'RowHeight', {34,'1x','1x'}, ...
    'ColumnWidth', {'1x','1x','1x'});
navigation = uigridlayout(qc,[1 4],'ColumnWidth',{36,36,'1x',80},'Padding',[0 0 0 0]);
navigation.Layout.Row = 1; navigation.Layout.Column = [1 2];
previousButton = uibutton(navigation,'Text','<','Tooltip','Previous sample (Left Arrow)', ...
    'Enable','off','Tag','PreviousSample','ButtonPushedFcn',@previousSample);
nextButton = uibutton(navigation,'Text','>','Tooltip','Next sample (Right Arrow)', ...
    'Enable','off','Tag','NextSample','ButtonPushedFcn',@nextSample);
sampleMenu = uidropdown(navigation, 'Items', {'No results loaded'}, 'Enable', 'off', ...
    'ValueChangedFcn', @showSample, 'Tag', 'SampleMenu');
sampleCounter = uilabel(navigation,'Text','0 / 0','HorizontalAlignment','center','Tag','SampleCounter');
qcActions = uigridlayout(qc,[1 2],'ColumnWidth',{'1x',115},'Padding',[0 0 0 0]);
qcActions.Layout.Row = 1; qcActions.Layout.Column = 3;
pdfButton = uibutton(qcActions, 'Text', 'Open QC PDF', 'Enable', 'off', 'ButtonPushedFcn', @openQC);
if ismac
    pdfButton.Text = 'Show QC PDF in Finder';
    pdfButton.Tooltip = 'Reveal the selected sample''s original QC PDF in Finder';
end
qcExportButton = uibutton(qcActions,'Text','Export PDFs...','Enable','off', ...
    'Tooltip','Regenerate QC PDFs for every eye in this loaded group without rerunning analysis', ...
    'Tag','ExportQCPDFs','ButtonPushedFcn',@exportQC);
imageAxes = uiaxes(qc); imageAxes.Layout.Row = 2; imageAxes.Layout.Column = 1;
angleAxes = uiaxes(qc); angleAxes.Layout.Row = 2; angleAxes.Layout.Column = 2;
medianAxes = uiaxes(qc); medianAxes.Layout.Row = 2; medianAxes.Layout.Column = 3;
profileAxes = uiaxes(qc); profileAxes.Layout.Row = 3; profileAxes.Layout.Column = [1 3];

results = uigridlayout(resultsTab, [4 2], 'RowHeight', {34,34,34,'1x'}, ...
    'ColumnWidth', {'1x',180});
resultLabel = uilabel(results, 'Text', 'No results loaded', 'WordWrap', 'on');
resultLabel.Layout.Row = 1; resultLabel.Layout.Column = [1 2];
folderButton = uibutton(results, 'Text', 'Open Output Folder', 'Enable', 'off', ...
    'ButtonPushedFcn', @(~,~) openFolder(fileparts(resultPath)));
folderButton.Layout.Row = 2; folderButton.Layout.Column = 2;
csvButton = uibutton(results, 'Text', 'Export Binned CSVs', 'Enable', 'off', 'ButtonPushedFcn', @exportCSV);
csvButton.Layout.Row = 3; csvButton.Layout.Column = 2;
binMenu = uidropdown(results, 'Items', {'1 degree','5 degrees','10 degrees'}, ...
    'Value', '5 degrees', 'ValueChangedFcn', @showBins);
binMenu.Layout.Row = 2; binMenu.Layout.Column = 1;
binTable = uitable(results, 'ColumnEditable', false);
binTable.Layout.Row = 4; binTable.Layout.Column = [1 2];
checkButton.Tooltip = 'Validate file matching, ROI geometry, image type and TIFF page; does not run analysis or assess anatomical quality.';
loadButton.Tooltip = 'Load a full group MAT into QC and Results without rerunning extraction.';
demoButton.Tooltip = 'Fill in example input paths and check them. Run Analysis or open saved results to populate QC.';
binMenu.Tooltip = 'Choose the saved angular summary table. Export Binned CSVs writes all three bin sizes.';
binTable.Tooltip = 'Numel counts pixels; STD is pixel spread; SE is for the mean; CI95L/U are pixel-based median limits. Click ? for definitions.';
controls = {roiField,imageField,outputField,groupField,channel,checkButton, ...
    runButton,loadButton,demoButton,folderButton,csvButton,pdfButton,qcExportButton,sampleMenu,binMenu};
toolboxNames = {'image_toolbox','curve_fitting_toolbox','statistics_toolbox'};
for n = 1:3
    if ~license('test', toolboxNames{n})
        environment.Text = 'Required toolbox license missing';
        environment.FontColor = [0.7 0.15 0.15];
    end
end
if ~isempty(options.Results.ResultsFile)
    readResults(options.Results.ResultsFile);
    tabs.SelectedTab = qcTab;
end

    function field = folderRow(row, text)
        label = uilabel(setup, 'Text', text); label.Layout.Row = row; label.Layout.Column = 1;
        field = uieditfield(setup, 'text'); field.Layout.Row = row; field.Layout.Column = 2;
        button = uibutton(setup, 'Text', 'Browse...', 'ButtonPushedFcn', @(~,~) browse(field));
        button.Layout.Row = row; button.Layout.Column = 3;
    end

    function browse(field)
        if busy, return; end
        selected = uigetdir(root, 'Select folder');
        if ~isequal(selected,0), field.Value = selected; end
    end

    function config = getConfig()
        config = struct('ROIdir',strtrim(roiField.Value), 'IMGdir',strtrim(imageField.Value), ...
            'Outputdir',strtrim(outputField.Value), 'GroupName',strtrim(groupField.Value), ...
            'Channel',channel.Value, 'Overwrite',false);
    end

    function useDemo(~,~)
        roiField.Value = fullfile(root,'JoVE_dataset','DMSO_4dpi','ROIs');
        imageField.Value = fullfile(root,'JoVE_dataset','DMSO_4dpi','TIFs 8-bit');
        outputField.Value = fullfile(root,'GUI_Output');
        groupField.Value = 'DMSO_4dpi'; channel.Value = 3;
        checkInputs([],[]);
    end

    function checkInputs(~,~)
        if busy, return; end
        try
            status.Text = 'Checking inputs...'; drawnow;
            manifest = rpegen_validate(getConfig());
            sampleTable.Data = manifest(:, {'Sample','Height','Width'});
            status.Text = sprintf('%d matched samples; inputs validated', height(manifest));
        catch problem
            sampleTable.Data = []; report(problem);
        end
    end

    function runAnalysis(~,~)
        if busy, return; end
        config = getConfig();
        try
            manifest = rpegen_validate(config);
            sampleTable.Data = manifest(:, {'Sample','Height','Width'});
            if isfile(fullfile(config.Outputdir, [config.GroupName '.mat']))
                choice = uiconfirm(fig, 'Replace this group''s MAT file and QC PDFs?', ...
                    'Existing Results', 'Options', {'Replace','Cancel'}, ...
                    'DefaultOption', 2, 'CancelOption', 2);
                if ~strcmp(choice,'Replace'), return; end
                config.Overwrite = true;
            end
            busy = true; setBusy(true);
            dialog = uiprogressdlg(fig, 'Title', 'RpEGEN', 'Message', 'Starting analysis', 'Value', 0);
            cleanup = onCleanup(@() finishRun(dialog));
            config.ProgressFcn = @progress;
            [newData, newPath] = rpegen_run(config);
            data = newData; resultPath = newPath;
            populate(); tabs.SelectedTab = qcTab;
            status.Text = sprintf('Complete: %d samples', numel(data.IMG));
        catch problem
            report(problem);
        end
        function progress(fraction, message)
            dialog.Value = fraction; dialog.Message = message;
            status.Text = sprintf('%s (%d%%)',message,round(fraction*100)); drawnow;
        end
        function finishRun(progressDialog)
            if isvalid(progressDialog), close(progressDialog); end
            busy = false; setBusy(false);
        end
    end

    function setBusy(value)
        comparison.SetBusy(value);
        for k = 1:numel(controls)
            controls{k}.Enable = matlab.lang.OnOffSwitchState(~value);
        end
        if ~value && isempty(data)
            resultControls = {folderButton,csvButton,pdfButton,qcExportButton,sampleMenu};
            for q = 1:numel(resultControls)
                resultControls{q}.Enable = 'off';
            end
        end
        updateNavigation();
    end

    function comparisonBusy(value)
        busy = value; setBusy(value);
    end

    function loadResults(~,~)
        if busy, return; end
        [name,folder] = uigetfile('*.mat', 'Select RpEGEN results', root);
        if isequal(name,0), return; end
        try
            readResults(fullfile(folder,name));
            tabs.SelectedTab = qcTab;
        catch problem
            report(problem);
        end
    end

    function readResults(file)
        loaded = load(file);
        names = fieldnames(loaded); candidate = [];
        required = {'IMG','IMGname','ROIxy','CLine','BIN_1_deg_all','BIN_5_deg_all','BIN_10_deg_all'};
        for k = 1:numel(names)
            entry = loaded.(names{k});
            if isstruct(entry) && isscalar(entry) && all(isfield(entry,required))
                candidate = entry; break;
            end
        end
        if isempty(candidate), error('RpEGEN:Results','This file does not contain RpEGEN results.'); end
        data = candidate; resultPath = char(file);
        populate(); status.Text = 'Results loaded';
    end

    function populate()
        sampleMenu.Items = data.IMGname(:)'; sampleMenu.Value = sampleMenu.Items{1};
        sampleMenu.Enable = 'on'; pdfButton.Enable = 'on'; qcExportButton.Enable = 'on';
        qcFolder = fileparts(resultPath);
        folderButton.Enable = 'on'; csvButton.Enable = 'on';
        [~,name] = fileparts(resultPath);
        resultLabel.Text = sprintf('%s | %d samples',name,numel(data.IMG));
        comparison.SetCurrent(resultPath);
        showSample([],[]); showBins([],[]);
    end

    function showHelp(~,~)
        topic = tabs.SelectedTab.Title;
        if tabs.SelectedTab==compareTab
            groupViews = findall(compareTab,'Tag','ComparisonViews');
            topic = groupViews.SelectedTab.Title;
        end
        [heading,body] = rpegen_help(topic);
        if isempty(helpWindow) || ~isvalid(helpWindow)
            helpWindow = uifigure('Name',heading,'Position',[180 140 660 570], ...
                'Tag','RpEGENHelp','Visible',fig.Visible);
            helpGrid = uigridlayout(helpWindow,[1 1]);
            uitextarea(helpGrid,'Editable','off','FontSize',14,'Tag','RpEGENHelpText');
        end
        helpWindow.Name = heading;
        helpText = findall(helpWindow,'Tag','RpEGENHelpText'); helpText.Value = body;
        if strcmp(fig.Visible,'on'), figure(helpWindow); end
    end

    function closeHelp()
        if ~isempty(helpWindow) && isvalid(helpWindow), delete(helpWindow); end
    end

    function showSample(~,~)
        if isempty(data), return; end
        updateNavigation();
        k = find(strcmp(data.IMGname,sampleMenu.Value),1);
        rpegen_draw_qc([imageAxes,angleAxes,medianAxes],profileAxes,data,k);
    end

    function previousSample(~,~)
        stepSample(-1);
    end

    function nextSample(~,~)
        stepSample(1);
    end

    function stepSample(step)
        if busy || isempty(data), return; end
        k = find(strcmp(data.IMGname,sampleMenu.Value),1);
        k = max(1,min(numel(data.IMGname),k+step));
        sampleMenu.Value = sampleMenu.Items{k}; showSample([],[]);
    end

    function updateNavigation()
        k = 0; count = 0;
        if ~isempty(data)
            count = numel(data.IMGname); k = find(strcmp(data.IMGname,sampleMenu.Value),1);
        end
        sampleCounter.Text = sprintf('%d / %d',k,count);
        previousButton.Enable = matlab.lang.OnOffSwitchState(~busy && k>1);
        nextButton.Enable = matlab.lang.OnOffSwitchState(~busy && k<count);
    end

    function sampleKey(~,event)
        if busy || ~isempty(event.Modifier), return; end
        % Leave native arrow-key editing and selection to input controls.
        active = fig.CurrentObject;
        if ismember(event.Key,{'leftarrow','rightarrow'}) && ~isempty(active) && ...
                (isprop(active,'Value') || isa(active,'matlab.ui.control.Table') || ...
                isa(active,'matlab.ui.container.TabGroup'))
            return;
        end
        if ismember(event.Key,{'leftarrow','pageup'})
            delta = -1;
        elseif ismember(event.Key,{'rightarrow','pagedown'})
            delta = 1;
        else
            return;
        end
        if tabs.SelectedTab==qcTab
            stepSample(delta);
        end
    end

    function showBins(~,~)
        if isempty(data), return; end
        fields = {'BIN_1_deg_all','BIN_5_deg_all','BIN_10_deg_all'};
        binTable.Data = data.(fields{find(strcmp(binMenu.Items,binMenu.Value),1)});
    end

    function openQC(~,~)
        target = fullfile(qcFolder, [sampleMenu.Value '_QC_plots.pdf']);
        if isfile(target)
            try
                if ismac
                    % Finder lets the user open the original with a direct file selection.
                    % A successful direct `open` launch can still be denied by Preview.
                    [code,message] = system(['/usr/bin/open -R ' shellQuote(target)]);
                    if code ~= 0
                        error('RpEGEN:OpenPDF', ...
                            'Could not reveal the QC PDF in Finder: %s\nFile:\n%s', ...
                            strtrim(message),target);
                    end
                else
                    open(target);
                end
            catch problem
                report(problem);
            end
        else
            uialert(fig,'No QC PDF exists for this sample.','QC PDF');
        end
    end

    function exportQC(~,~)
        if busy || isempty(data), return; end
        folder = uigetdir(qcFolder,'Export QC PDFs for this group');
        if isequal(folder,0), return; end
        targets = fullfile(folder,strcat(string(data.IMGname),'_QC_plots.pdf'));
        if any(isfile(targets))
            choice = uiconfirm(fig,'Replace existing QC PDFs in this folder? Source images and MAT results will not change.', ...
                'Existing QC PDFs','Options',{'Replace','Cancel'},'DefaultOption',2,'CancelOption',2);
            if ~strcmp(choice,'Replace'), return; end
        end
        dialog = uiprogressdlg(fig,'Title','Export QC PDFs','Message','Starting export','Value',0);
        busy = true; setBusy(true);
        cleanup = onCleanup(@finishExport);
        try
            files = rpegen_export_qc(data,folder,@exportProgress);
            qcFolder = folder;
            status.Text = sprintf('Exported %d QC PDFs to %s',numel(files),folder);
        catch problem
            report(problem);
        end
        function exportProgress(fraction,message)
            dialog.Value = fraction; dialog.Message = message; drawnow;
        end
        function finishExport()
            if isvalid(dialog), close(dialog); end
            busy = false; setBusy(false);
        end
    end

    function exportCSV(~,~)
        folder = uigetdir(fileparts(resultPath),'Select CSV export folder');
        if isequal(folder,0), return; end
        [~,name] = fileparts(resultPath);
        fields = {'BIN_1_deg_all','BIN_5_deg_all','BIN_10_deg_all'};
        targets = cellfun(@(field) fullfile(folder,[name '_' field '.csv']),fields,'UniformOutput',false);
        if any(cellfun(@isfile,targets))
            choice = uiconfirm(fig,'Replace existing binned CSV files?','Existing CSVs', ...
                'Options',{'Replace','Cancel'},'DefaultOption',2,'CancelOption',2);
            if ~strcmp(choice,'Replace'), return; end
        end
        try
            for k = 1:numel(fields), writetable(data.(fields{k}),targets{k}); end
            status.Text = 'Binned CSV files exported';
        catch problem
            report(problem);
        end
    end

    function report(problem)
        status.Text = 'Action failed';
        fprintf(2, '%s\n', getReport(problem, 'extended', 'hyperlinks', 'off'));
        uialert(fig,problem.message,'RpEGEN');
    end

    function closeApp(~,~)
        if busy
            uialert(fig,'A calculation is running. Finish it or cancel the permutation test before closing.','Calculation Running');
        else
            delete(fig);
        end
    end
end

function openFolder(folder)
if ispc
    winopen(folder);
elseif ismac
    system(['open ' shellQuote(folder)]);
else
    system(['xdg-open ' shellQuote(folder)]);
end
end

function quoted = shellQuote(value)
quoted = ['''' strrep(value, '''', '''"''"''') ''''];
end
