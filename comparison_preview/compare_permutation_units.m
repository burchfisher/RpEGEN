function compare_permutation_units
% Exploratory comparison only: no edits to permsim, extraction, or source data.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
names = {'DMSO_4dpi','IWR1_4dpi'}; datasets = cell(1,2); sizes = zeros(1,2);
for g = 1:2
    loaded = load(fullfile(root,'JoVE_dataset','Output',[names{g} '.mat']),names{g});
    datasets{g} = loaded.(names{g}); sizes(g) = numel(datasets{g}.IMGname);
end
iterations = 5000; seed = 20261013;
previous = rng; cleanupRng = onCleanup(@() rng(previous)); rng(seed,'twister');
verifyPixelSampler();
count = sum(sizes); membership = zeros(count,iterations);
for draw = 1:iterations
    membership(randperm(count,sizes(1)),draw) = 1;
end
assert(all(sum(membership,1) == sizes(1)));
angles = (2.5:5:177.5)'; medians = zeros(36,2);
pixelExceedances = zeros(36,1); eyeExceedances = pixelExceedances;
for b = 1:36
    histograms = zeros(256,count); offset = 0;
    for g = 1:2
        data = datasets{g}; raw = data.RAW_data;
        [known,ids] = ismember(string(raw.IMG_name),string(data.IMGname(:)));
        assert(all(known));
        inBin = raw.Angle >= (b-1)*5 & raw.Angle < b*5;
        values = double(raw.Raw_Data(inBin));
        assert(all(values >= 0 & values <= 255 & values == fix(values)));
        histograms(:,offset+(1:sizes(g))) = accumarray([values+1,ids(inBin)],1,[256 sizes(g)]);
        medians(b,g) = median(values);
        assert(medians(b,g) == data.BIN_5_deg_all.Median(b));
        offset = offset+sizes(g);
    end
    assert(all(sum(histograms,1) > 0));
    difference = medians(b,1)-medians(b,2);
    totalHistogram = sum(histograms,2);
    groupAHistogram = sum(histograms(:,1:sizes(1)),2);
    % Literal pixel shuffles reproduce permsim, using memory-bounded batches.
    pixelExceedances(b) = pixelShuffleExceedances(totalHistogram,sum(groupAHistogram),iterations,difference);
    % Reuse these complete-eye assignments at every angular bin.
    eyeA = histograms*membership;
    eyeB = totalHistogram-eyeA;
    eyeDifferences = histogramMedian(eyeA)-histogramMedian(eyeB);
    eyeExceedances(b) = nnz(abs(eyeDifferences) >= abs(difference));
    for draw = 1:3
        literalA = repelem((0:255)',eyeA(:,draw));
        literalB = repelem((0:255)',eyeB(:,draw));
        assert(median(literalA)-median(literalB) == eyeDifferences(draw));
    end
    fprintf('Bin %2d/36 (%3d-%3d): difference=%6.1f; pixel p=%.5f; eye p=%.5f\n', ...
        b,(b-1)*5,b*5,difference,(pixelExceedances(b)+1)/(iterations+1), ...
        (eyeExceedances(b)+1)/(iterations+1));
end
pPixel = (pixelExceedances+1)/(iterations+1);
pEye = (eyeExceedances+1)/(iterations+1);
legacyUncorrectedMonteCarloP = pixelExceedances/iterations;
holmPixel = holmAdjustment(pPixel); holmEye = holmAdjustment(pEye);
assert(all(holmPixel >= pPixel) && all(holmEye >= pEye));
fprintf('\nUnadjusted p <= .05: pixel %d/36, eye %d/36\n',nnz(pPixel<=.05),nnz(pEye<=.05));
fprintf('Holm adjusted across 36 bins: pixel %d/36, eye %d/36\n',nnz(holmPixel<=.05),nnz(holmEye<=.05));
printRanges('Eye shuffle, unadjusted',pEye<=.05);
printRanges('Eye shuffle, Holm-adjusted',holmEye<=.05);
fprintf('Eye p-value range: %.5f..%.5f; minimum Holm eye p: %.5f\n',min(pEye),max(pEye),min(holmEye));
fprintf('Eye-shuffle bins within .01 of raw .05 threshold: '); fprintf('%.1f ',angles(abs(pEye-.05)<.01)); fprintf('\n');

fig = figure('Visible','off','Color','w','Position',[100 100 1220 1050]);
cleanup = onCleanup(@() close(fig));
tiles = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(tiles); hold(ax,'on');
plot(ax,angles,medians(:,1),'Color',[0 .40 .70],'LineWidth',2);
plot(ax,angles,medians(:,2),'Color',[.7 .12 .35],'LineWidth',2);
xlim(ax,[0 180]); ylim(ax,[0 255]); grid(ax,'on'); ylabel(ax,'Median intensity (0-255)');
title(ax,'Same observed pooled-pixel medians; only the permutation unit changes');
legend(ax,{sprintf('DMSO_4dpi (n = %d)',sizes(1)),sprintf('IWR1_4dpi (n = %d)',sizes(2))}, ...
    'Interpreter','none','Location','eastoutside'); ax.FontSize = 11;
drawP(nexttile(tiles),angles,pPixel,pEye,iterations,'Per-angle p-values (unadjusted)');
drawP(nexttile(tiles),angles,holmPixel,holmEye,iterations,'P-values adjusted for all 36 angle bins (Holm)');
title(tiles,'DMSO_4dpi vs IWR1_4dpi | 5-degree bins | 5,000 permutations','Interpreter','none','FontSize',16);
xlabel(tiles,{'Exploratory whole-eye test assumes exchangeable, independent eyes; wells/batches are not modeled.', ...
    'Two-sided median difference; Monte Carlo p = (exceedances + 1)/(5,000 + 1). Lowest plotted p is the simulation limit, not zero.'}, ...
    'FontSize',10);
destination = fullfile(root,'comparison_preview','example_permutation_comparison.png');
exportgraphics(fig,destination,'Resolution',170);
save(fullfile(root,'comparison_preview','example_permutation_comparison.mat'), ...
    'names','sizes','angles','iterations','seed','medians','pPixel','pEye','holmPixel','holmEye', ...
    'pixelExceedances','eyeExceedances','legacyUncorrectedMonteCarloP','membership');
fprintf('Saved: %s\n',destination);
end

function exceedances = pixelShuffleExceedances(counts,n,iterations,observed)
population = repelem((0:numel(counts)-1)',counts);
total = numel(population); exceedances = 0;
for first = 1:128:iterations
    batch = min(128,iterations-first+1); indices = zeros(total,batch);
    for draw = 1:batch, indices(:,draw) = randperm(total)'; end
    randomized = population(indices);
    differences = median(randomized(1:n,:),1)-median(randomized(n+1:end,:),1);
    exceedances = exceedances+nnz(abs(differences) >= abs(observed));
end
end

function values = histogramMedian(counts)
cumulative = cumsum(counts,1); total = cumulative(end,:);
values = (sum(cumulative < ceil(total/2),1)+sum(cumulative < floor(total/2)+1,1))/2;
end

function adjusted = holmAdjustment(p)
[ordered,index] = sort(p); m = numel(p);
adjusted = zeros(size(p)); adjusted(index) = min(1,cummax((m:-1:1)'.*ordered));
end

function verifyPixelSampler()
% The exact same RNG sequence must reproduce the original implementation.
previous = rng; cleanup = onCleanup(@() rng(previous));
population = [0;0;1;2;3;3;4]; counts = accumarray(population+1,1,[5 1]);
for n = [2 3 4]
    observed = median(population(1:n))-median(population(n+1:end));
    rng(17,'twister'); reference = permsim(population(1:n),population(n+1:end),257);
    rng(17,'twister'); actual = pixelShuffleExceedances(counts,n,257,observed)/257;
    assert(abs(reference-actual) < 1e-12,'Chunked pixel shuffles differ from permsim.');
end
assert(isequal(holmAdjustment([.01;.02;.2]),[.03;.04;.2]));
end

function drawP(ax,angles,pPixel,pEye,iterations,heading)
semilogy(ax,angles,pPixel,'Color',[.25 .25 .25],'LineWidth',1.8); hold(ax,'on');
semilogy(ax,angles,pEye,'Color',[0 .48 .36],'LineWidth',2);
yline(ax,.05,'--','0.05','Color',[.72 .16 .16],'LabelHorizontalAlignment','left');
xlim(ax,[0 180]); ylim(ax,[1/(iterations+1) 1]); grid(ax,'on');
ax.YTick = [1/(iterations+1) .001 .01 .05 .1 1];
ax.YTickLabel = {'0.0002','0.001','0.01','0.05','0.1','1'};
ylabel(ax,'p-value (log scale)'); xlabel(ax,'Angular position (degrees)');
title(ax,heading);
legend(ax,{'Pixel shuffle (legacy-style)','Whole-eye shuffle'},'Location','eastoutside');
ax.FontSize = 11;
end

function printRanges(label,significant)
edges = diff([false;significant;false]); starts = find(edges==1); stops = find(edges==-1)-1;
fprintf('%s intervals: ',label);
for k = 1:numel(starts), fprintf('[%d, %d) ',(starts(k)-1)*5,stops(k)*5); end
if isempty(starts), fprintf('none'); end
fprintf('\n');
end
