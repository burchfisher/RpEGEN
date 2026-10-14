function test_rpegen_extraction
% Re-extract bundled inputs into a fresh folder and compare with saved results.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
output = tempname; mkdir(output);
config = struct('ROIdir',fullfile(root,'JoVE_dataset','DMSO_4dpi','ROIs'), ...
    'IMGdir',fullfile(root,'JoVE_dataset','DMSO_4dpi','TIFs 8-bit'), ...
    'Outputdir',output,'GroupName','DMSO_4dpi','Channel',3,'Overwrite',false);
beforeFolder = pwd; beforePath = path;
[data,file] = rpegen_run(config);
assert(strcmp(pwd,beforeFolder) && strcmp(path,beforePath));
saved = load(fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'));
assert(numel(data.IMGname)==11 && isequal(data.IMGname,saved.DMSO_4dpi.IMGname));
for field = {'RAW_data','BIN_1_deg_all','BIN_5_deg_all','BIN_10_deg_all'}
    assert(isequaln(data.(field{1}),saved.DMSO_4dpi.(field{1})), ...
        'Extraction differs from bundled reference: %s',field{1});
end
assert(isfile(file) && isfile(fullfile(output,'DMSO_4dpi_run.mat')));
assert(numel(dir(fullfile(output,'*_QC_plots.pdf')))==11);
record = load(fullfile(output,'DMSO_4dpi_run.mat'));
assert(height(record.manifest)==11 && record.settings.Channel==3);
try
    rpegen_run(config);
    error('Test:ExpectedError','Existing output was overwritten without permission.');
catch problem
    assert(strcmp(problem.identifier,'RpEGEN:Exists'));
end
fprintf('Fresh extraction passed: 11 eyes, exact raw/table parity, QC PDFs, run record, and overwrite protection. Outputs: %s\n',output);
end
