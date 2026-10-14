function compare_pixel_counts
% Descriptive counts and sensitivity, with no changes to analysis or input data.
root = fileparts(fileparts(mfilename('fullpath')));
source = fullfile(root,'GUI_Output','DMSO_4dpi.mat');
if ~isfile(source)
    source = fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat');
end
loaded = load(source,'DMSO_4dpi'); data = loaded.DMSO_4dpi;
raw = data.RAW_data; eyeNames = string(data.IMGname(:));
[known,eyeIndex] = ismember(string(raw.IMG_name),eyeNames);
assert(all(known) && numel(unique(eyeNames)) == numel(eyeNames));
eyeCount = numel(eyeNames); angles = (2.5:5:177.5)';
counts = zeros(36,eyeCount); pooledMedian = nan(36,1);
leaveOneOut = nan(36,eyeCount);
for b = 1:36
    inBin = raw.Angle >= (b-1)*5 & raw.Angle < b*5;
    pooledMedian(b) = median(double(raw.Raw_Data(inBin)));
    for e = 1:eyeCount
        counts(b,e) = nnz(inBin & eyeIndex == e);
        leaveOneOut(b,e) = median(double(raw.Raw_Data(inBin & eyeIndex ~= e)));
    end
end
assert(isequal(sum(counts,2),data.BIN_5_deg_all.Numel));
assert(isequaln(pooledMedian,data.BIN_5_deg_all.Median));
meanCount = mean(counts,2);
sdCount = std(counts,0,2); % Sample SD across eyes, with N-1 denominator.
shares = 100*counts./sum(counts,2);
assert(all(abs(sum(shares,2)-100) < 1e-10));
[maxShare,largestEye] = max(shares,[],2);
[peakShare,peakBin] = max(maxShare);
absoluteChange = abs(leaveOneOut-pooledMedian);
[maxChange,influentialEye] = max(absoluteChange,[],2);
[peakChange,changeBin] = max(maxChange);
totalShares = 100*sum(counts,1)/sum(counts(:));
[maxTotalShare,totalEye] = max(totalShares);
shortNames = eyeNames;
hasSuffix = contains(shortNames,' central');
shortNames(hasSuffix) = extractBefore(shortNames(hasSuffix),' central');
fprintf('Source: %s\n',source);
fprintf('%d eyes; equal contribution would be %.2f%% in each bin.\n',eyeCount,100/eyeCount);
fprintf('Mean count across eyes ranges %.1f..%.1f; sample SD ranges %.1f..%.1f\n', ...
    min(meanCount),max(meanCount),min(sdCount),max(sdCount));
fprintf('Largest contributor per bin: median %.2f%%, range %.2f..%.2f%%\n', ...
    median(maxShare),min(maxShare),max(maxShare));
fprintf('Peak individual share %.2f%% at %d-%d degrees: %s (%d of %d pixels)\n', ...
    peakShare,(peakBin-1)*5,peakBin*5,shortNames(largestEye(peakBin)), ...
    counts(peakBin,largestEye(peakBin)),sum(counts(peakBin,:)));
fprintf('Largest overall pixel share: %.2f%% from %s\n',maxTotalShare,shortNames(totalEye));
fprintf('At 40-45 degrees: mean %.2f, SD %.2f, min %d, max %d, largest share %.2f%%\n', ...
    meanCount(9),sdCount(9),min(counts(9,:)),max(counts(9,:)),maxShare(9));
fprintf('Largest leave-one-eye-out median change %.2f at %d-%d degrees, omitting %s (pixel share %.2f%%)\n', ...
    peakChange,(changeBin-1)*5,changeBin*5,shortNames(influentialEye(changeBin)), ...
    shares(changeBin,influentialEye(changeBin)));
fprintf('Median of maximum leave-one-eye-out changes across bins: %.2f\n',median(maxChange));
fprintf('Highest contributor counts (number of angular bins):\n');
for e = 1:eyeCount
    fprintf('  %s: %d bins; total share %.2f%%\n',shortNames(e),nnz(largestEye==e),totalShares(e));
end

fig = figure('Visible','off','Color','w','Position',[100 100 1450 1050]);
cleanup = onCleanup(@() close(fig));
tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tiles,sprintf('DMSO_4dpi | %d eyes | 5-degree bins',eyeCount), ...
    'Interpreter','none','FontSize',17,'FontWeight','bold');
ax = nexttile(tiles); hold(ax,'on');
spread = fill(ax,[angles;flipud(angles)], ...
    [meanCount-sdCount;flipud(meanCount+sdCount)],[.08 .55 .40], ...
    'FaceAlpha',.16,'EdgeColor','none');
individuals = plot(ax,angles,counts,'Color',[.72 .72 .72],'LineWidth',.8);
average = plot(ax,angles,meanCount,'Color',[0 .40 .28],'LineWidth',2.5);
title(ax,'Pixel count per eye: mean and sample SD');
ylabel(ax,'Pixels per eye per bin'); xlabel(ax,'Angular position (degrees)');
xlim(ax,[0 180]); ylim(ax,[0 ceil(max(counts(:))/100)*100+100]); grid(ax,'on');
legend(ax,[average,spread,individuals(1)],{'Mean','Mean +/- 1 SD','Individual eyes'}, ...
    'Location','southoutside','NumColumns',3);
ax.FontSize = 11;

shareAx = nexttile(tiles); hold(shareAx,'on');
plot(shareAx,angles,maxShare,'Color',[.80 .29 .06],'LineWidth',2.4);
yline(shareAx,100/eyeCount,'--','Equal share (1/11)','Color',[.3 .3 .3], ...
    'LabelHorizontalAlignment','left');
plot(shareAx,angles(peakBin),peakShare,'o','MarkerFaceColor',[.80 .29 .06], ...
    'Color',[.80 .29 .06],'MarkerSize',6);
title(shareAx,'Largest individual contribution at each angle');
ylabel(shareAx,'Share of pooled pixels (%)'); xlabel(shareAx,'Angular position (degrees)');
xlim(shareAx,[0 180]); ylim(shareAx,[0 max(30,ceil(peakShare/5)*5+5)]); grid(shareAx,'on');
shareAx.FontSize = 11;

heatAx = nexttile(tiles);
imagesc(heatAx,angles,1:eyeCount,shares');
heatAx.YTick = 1:eyeCount; heatAx.YTickLabel = shortNames;
heatAx.TickLabelInterpreter = 'none';
xlim(heatAx,[0 180]); clim(heatAx,[0 ceil(peakShare/5)*5]); colormap(heatAx,parula(256));
bar = colorbar(heatAx); bar.Label.String = 'Share of pooled pixels (%)';
title(heatAx,'Which eye contributes the pixels?');
xlabel(heatAx,'Angular position (degrees)'); heatAx.FontSize = 10;

influenceAx = nexttile(tiles); hold(influenceAx,'on');
plot(influenceAx,angles,maxChange,'Color',[.67 .12 .38],'LineWidth',2.4);
title(influenceAx,'Sensitivity: remove one eye at a time');
subtitle(influenceAx,'Largest absolute change in the pooled median at each bin');
xlabel(influenceAx,'Angular position (degrees)'); ylabel(influenceAx,'Median change (intensity units)');
xlim(influenceAx,[0 180]); ylim(influenceAx,[0 max(5,ceil(peakChange/5)*5+5)]); grid(influenceAx,'on');
influenceAx.FontSize = 11;
xlabel(tiles,{'Counts include every eye in every bin; shading is sample SD, not a confidence interval.', ...
    'The largest contributor can be a different eye at each angle. Leave-one-out change is a descriptive sensitivity check.'}, ...
    'FontSize',11);
destination = fullfile(root,'comparison_preview','DMSO_4dpi_pixel_counts.png');
exportgraphics(fig,destination,'Resolution',160);
save(fullfile(root,'comparison_preview','DMSO_4dpi_pixel_counts.mat'), ...
    'source','angles','eyeNames','counts','meanCount','sdCount','shares','maxShare', ...
    'largestEye','pooledMedian','leaveOneOut','maxChange','influentialEye');
fprintf('Preview: %s\n',destination);
end
