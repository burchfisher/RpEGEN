function test_rpegen_qc_and_batch_exports
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
source = fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat');
loaded = load(source); data = loaded.DMSO_4dpi;
for field = {'IMG','IMGname','ROIxy','ROIBW','CLine'}
    data.(field{1}) = data.(field{1})(1:2);
end
original = data;
folder = fullfile(tempdir,'rpegen_qc_batch'); if ~isfolder(folder), mkdir(folder); end
progressValues = [];
files = rpegen_export_qc(data,folder,@recordProgress);
assert(numel(files)==2 && all(isfile(files)));
assert(progressValues(1)==0 && progressValues(end)==1 && all(diff(progressValues)>=0));
assert(isempty(findall(groot,'Tag','RpEGENExport')));
assert(isequaln(data,original));

fig = RpEGEN_GUI('Visible','off','ResultsFile',source);
cleanup = onCleanup(@() delete(fig));
axesHandles = [findall(fig,'Tag','QCImageAxes'),findall(fig,'Tag','QCAngleAxes'),findall(fig,'Tag','QCMedianAxes')];
assert(numel(axesHandles)==3);
for j = 1:3
    ax = axesHandles(j); image = findall(ax,'Type','image');
    assert(isscalar(image) && isequal(image.CData,repmat(double(data.IMG{1})/255,1,1,3)));
    assert(strcmp(ax.YDir,'reverse'));
    assert(isequal(ax.XLim,[.5 size(data.IMG{1},2)+.5]) && isequal(ax.YLim,[.5 size(data.IMG{1},1)+.5]));
end
assert(numel(findall(fig,'Type','colorbar'))==3);
points = findall(fig,'Tag','QCMedianPoints');
assert(isequal(points.CData(:),data.CLine{1}.median(:)) && points.SizeData==10);
assert(isequal(points.YData(:),data.CLine{1}.y(:)));
assert(strcmp(findall(fig,'Tag','ExportQCPDFs').Enable,'on'));
for viewSize = [1180 780;1000 700]'
    fig.Position(3:4) = viewSize'; fig.Visible = 'on'; drawnow; pause(.5); drawnow;
    positions = zeros(3,4);
    for j = 1:3, axesHandles(j).Units = 'pixels'; positions(j,:) = tightPosition(axesHandles(j)); end
    assert(max(positions(:,3))-min(positions(:,3))<2,'QC image widths differ.');
    assert(max(positions(:,4))-min(positions(:,4))<2,'QC image heights differ.');
end
exportapp(fig,fullfile(folder,'qc_gui.png'));
clear cleanup;

[~,help] = rpegen_help('Median Profiles');
assert(any(contains(help,'A, A, C')) && any(contains(help,'That is the pooled median')));
fprintf('QC overlays, equal panel sizes, saved-data exports, and CI help passed: %s\n',folder);
    function recordProgress(value,~)
        progressValues(end+1) = value;
    end
end
