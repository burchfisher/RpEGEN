function compare_median_definitions
% Descriptive preview only. Does not modify analysis code or source MAT files.
root = fileparts(fileparts(mfilename('fullpath')));
source = fullfile(root,'GUI_Output','DMSO_4dpi.mat');
if ~isfile(source)
    source = fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat');
end
loaded = load(source,'DMSO_4dpi'); data = loaded.DMSO_4dpi;
raw = data.RAW_data;
eyeNames = string(data.IMGname(:));
assert(numel(unique(eyeNames)) == numel(eyeNames));
[known,eyeIndex] = ismember(string(raw.IMG_name),eyeNames);
assert(all(known),'Raw measurements contain unknown eye identifiers.');
angles = (2.5:5:177.5)';
eyeMedians = nan(36,numel(eyeNames));
pixelCounts = zeros(size(eyeMedians));
pixelHistograms = zeros(256,numel(eyeNames),36);
pooledMedian = nan(36,1);
for b = 1:36
    % Match the published half-open bins, including exclusion of angle 180.
    inBin = raw.Angle >= (b-1)*5 & raw.Angle < b*5;
    pooledMedian(b) = median(double(raw.Raw_Data(inBin)));
    for e = 1:numel(eyeNames)
        values = double(raw.Raw_Data(inBin & eyeIndex == e));
        pixelCounts(b,e) = numel(values);
        if ~isempty(values)
            assert(all(values >= 0 & values <= 255 & values == fix(values)));
            eyeMedians(b,e) = median(values);
            pixelHistograms(:,e,b) = accumarray(values+1,1,[256 1]);
        end
    end
end
assert(isequaln(pooledMedian,data.BIN_5_deg_all.Median));
assert(all(pixelCounts(:) > 0),'Some eyes do not contribute to all bins.');
equalEyeMedian = median(eyeMedians,2);
difference = equalEyeMedian-pooledMedian;
[maxAbsoluteDifference,maxBin] = max(abs(difference));
fprintf('Source: %s\n',source);
fprintf('%d eyes; all eyes contribute to every 5-degree bin.\n',numel(eyeNames));
fprintf('Median absolute difference: %.2f intensity units\n',median(abs(difference)));
fprintf('Mean absolute difference: %.2f intensity units\n',mean(abs(difference)));
fprintf('Maximum absolute difference: %.2f at %.1f degrees (bin %d-%d)\n', ...
    maxAbsoluteDifference,angles(maxBin),(maxBin-1)*5,maxBin*5);
fprintf('At this bin: pooled=%.2f, median of eye medians=%.2f, difference=%+.2f\n', ...
    pooledMedian(maxBin),equalEyeMedian(maxBin),difference(maxBin));
fprintf('Pixels per eye in this bin: '); fprintf('%d ',pixelCounts(maxBin,:)); fprintf('\n');
fprintf('Individual-eye medians in this bin: '); fprintf('%.1f ',eyeMedians(maxBin,:)); fprintf('\n');

% Resample the same whole eyes for every angular position and both summaries.
% Histograms reproduce the exact median of duplicated 8-bit pixel values.
previousRng = rng; rngCleanup = onCleanup(@() rng(previousRng));
rng(20261012,'twister');
iterations = 5000; eyeCount = numel(eyeNames);
[~,sampledEyes] = bootstrp(iterations,@mean,(1:eyeCount)');
drawNumber = repmat(1:iterations,eyeCount,1);
weights = accumarray([sampledEyes(:),drawNumber(:)],1,[eyeCount iterations]);
pooledBoot = nan(36,iterations); equalEyeBoot = pooledBoot;
for b = 1:36
    pooledBoot(b,:) = histogramMedian(pixelHistograms(:,:,b)*weights);
    equalEyeBoot(b,:) = median(reshape(eyeMedians(b,sampledEyes(:)),eyeCount,iterations),1);
    assert(histogramMedian(sum(pixelHistograms(:,:,b),2)) == pooledMedian(b));
    % Cross-check the optimized bootstrap against literal duplicated pixels.
    inBin = raw.Angle >= (b-1)*5 & raw.Angle < b*5;
    for draw = 1:3
        values = cell(eyeCount,1);
        for e = 1:eyeCount
            values{e} = double(raw.Raw_Data(inBin & eyeIndex == sampledEyes(e,draw)));
        end
        assert(median(vertcat(values{:})) == pooledBoot(b,draw));
    end
end
pooledEyeCI = prctile(pooledBoot,[2.5 97.5],2);
equalEyeCI = prctile(equalEyeBoot,[2.5 97.5],2);
legacyCI = [data.BIN_5_deg_all.CI95L,data.BIN_5_deg_all.CI95U];
assert(all(pooledEyeCI(:,1) <= pooledEyeCI(:,2)));
assert(all(equalEyeCI(:,1) <= equalEyeCI(:,2)));
fprintf('Median CI widths: legacy pixels %.2f; pooled/eye-bootstrap %.2f; equal-eye/eye-bootstrap %.2f\n', ...
    median(diff(legacyCI,1,2)),median(diff(pooledEyeCI,1,2)),median(diff(equalEyeCI,1,2)));
fprintf('At %.1f degrees: legacy CI [%.2f, %.2f]; pooled eye-bootstrap CI [%.2f, %.2f]; equal-eye CI [%.2f, %.2f]\n', ...
    angles(maxBin),legacyCI(maxBin,:),pooledEyeCI(maxBin,:),equalEyeCI(maxBin,:));

fig = figure('Visible','off','Color','w','Position',[100 100 1250 1000]);
cleanup = onCleanup(@() close(fig));
tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tiles,[1 2]); hold(ax,'on');
band(ax,angles,legacyCI,[0 .40 .70],.25);
band(ax,angles,equalEyeCI,[.75 .12 .34],.18);
pooled = plot(ax,angles,pooledMedian,'Color',[0 .40 .70],'LineWidth',2.6);
equal = plot(ax,angles,equalEyeMedian,'Color',[.75 .12 .34],'LineWidth',2.6);
xlim(ax,[0 180]); ylim(ax,[0 255]); grid(ax,'on');
ylabel(ax,'Median pixel intensity (0-255)');
title(ax,sprintf('DMSO_4dpi | %d eyes | 5-degree bins',numel(eyeNames)), ...
    'Interpreter','none','FontSize',16);
subtitle(ax,'Current pooled-pixel median vs median of per-eye medians');
legend(ax,[pooled,equal], ...
    {'Pooled pixels + nominal 95% pixel-based median CI','Equal-eye summary + 95% eye-bootstrap CI'}, ...
    'Location','southoutside','NumColumns',2);
ax.FontSize = 12;
controlAx = nexttile(tiles); hold(controlAx,'on');
eyeBand = band(controlAx,angles,pooledEyeCI,[0 .40 .70],.18);
legacyBand = band(controlAx,angles,legacyCI,[0 .40 .70],.60);
plot(controlAx,angles,pooledMedian,'Color',[0 .40 .70],'LineWidth',2);
xlim(controlAx,[0 180]); ylim(controlAx,[0 255]); grid(controlAx,'on');
title(controlAx,'Same pooled median, different CI method');
ylabel(controlAx,'Median pixel intensity (0-255)');
xlabel(controlAx,'Angular position (degrees)');
legend(controlAx,[legacyBand,eyeBand],{'Pixel-based median CI','Whole-eye bootstrap CI'}, ...
    'Location','southoutside','NumColumns',2);
controlAx.FontSize = 11;
deltaAx = nexttile(tiles); hold(deltaAx,'on');
yline(deltaAx,0,'Color',[.4 .4 .4],'LineWidth',1);
plot(deltaAx,angles,difference,'Color',[.75 .12 .34],'LineWidth',2);
plot(deltaAx,angles(maxBin),difference(maxBin),'o','Color',[.75 .12 .34], ...
    'MarkerFaceColor',[.75 .12 .34],'MarkerSize',7);
limit = max(5,ceil(maxAbsoluteDifference/5)*5+5);
xlim(deltaAx,[0 180]); ylim(deltaAx,[-limit limit]); grid(deltaAx,'on');
xlabel(deltaAx,'Angular position (degrees)');
ylabel(deltaAx,'Difference (intensity units)');
title(deltaAx,'Median of eye medians minus pooled-pixel median');
deltaAx.FontSize = 11;
xlabel(tiles,{'5,000 whole-eye resamples; pointwise 95% percentile bootstrap intervals.', ...
    'Illustrative: assumes independent eyes; treatment wells and batches are not modeled. Dorsal 0, ventral 180.'}, ...
    'FontSize',11);
destination = fullfile(root,'comparison_preview','DMSO_4dpi_median_comparison.png');
exportgraphics(fig,destination,'Resolution',170);
save(fullfile(root,'comparison_preview','DMSO_4dpi_median_comparison.mat'), ...
    'source','angles','eyeNames','eyeMedians','pixelCounts','pooledMedian', ...
    'equalEyeMedian','difference','legacyCI','pooledEyeCI','equalEyeCI','iterations','sampledEyes');
fprintf('Preview: %s\n',destination);
end

function values = histogramMedian(counts)
cumulative = cumsum(counts,1); total = cumulative(end,:);
lower = sum(cumulative < ceil(total/2),1);
upper = sum(cumulative < floor(total/2)+1,1);
values = (lower+upper)/2;
end

function h = band(ax,x,limits,color,alpha)
h = fill(ax,[x;flipud(x)],[limits(:,1);flipud(limits(:,2))],color, ...
    'FaceAlpha',alpha,'EdgeColor','none');
end
