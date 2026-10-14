function [data, outputFile] = rpegen_run(config)
% Run the original pipeline in an isolated function workspace.
root = fileparts(mfilename('fullpath'));
oldFolder = pwd;
oldPath = path;
cleanup = onCleanup(@() restoreEnvironment(oldFolder, oldPath));
addpath(root, fullfile(root, 'ReadImageJROI-master'));
if ~isfield(config, 'ProgressFcn')
    config.ProgressFcn = @(fraction, message) [];
end
for toolbox = {'image_toolbox', 'curve_fitting_toolbox', 'statistics_toolbox'}
    if ~license('test', toolbox{1})
        error('RpEGEN:Toolbox', 'Required toolbox license unavailable: %s', toolbox{1});
    end
end
manifest = rpegen_validate(config);
if isempty(config.Outputdir)
    error('RpEGEN:Output', 'Select an output folder.');
end
if ~isfolder(config.Outputdir)
    [ok, message] = mkdir(config.Outputdir);
    if ~ok, error('RpEGEN:Output', '%s', message); end
end
outputFile = fullfile(config.Outputdir, [config.GroupName '.mat']);
if isfile(outputFile) && (~isfield(config, 'Overwrite') || ~config.Overwrite)
    error('RpEGEN:Exists', 'Results already exist: %s. Choose another group name or allow replacement.', outputFile);
end
rpegenConfig = config;
data = struct;
rpegenConfig.ProgressFcn(0.05, 'Loading images and ROIs');
logText = evalc('run(fullfile(root, ''RpEGEN.m''))');
if ~exist('data', 'var') || ~isfield(data, 'BIN_10_deg_all')
    error('RpEGEN:Incomplete', 'Analysis did not complete.\n%s', logText);
end
settings = rmfield(config, 'ProgressFcn');
save(fullfile(config.Outputdir, [config.GroupName '_run.mat']), 'settings', 'manifest', 'logText');
config.ProgressFcn(1, 'Analysis complete');
end

function restoreEnvironment(oldFolder, oldPath)
cd(oldFolder);
path(oldPath);
close(findall(groot, 'Type', 'figure', 'Tag', 'RpEGENExport'));
end
