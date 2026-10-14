function audit_permutation_precision
% Audit the existing pixel-level test, without changing extraction or MATs.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder = fullfile(root,'comparison_preview');
left = load(fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'));
right = load(fullfile(root,'JoVE_dataset','Output','IWR1_4dpi.mat'));
a = left.DMSO_4dpi.RAW_data; b = right.IWR1_4dpi.RAW_data;
previous = rng; cleanup = onCleanup(@() rng(previous));
plan = rpegen_permutation_plan;
angle = (0.5:1:179.5)'; p5000 = zeros(180,1); mc = zeros(180,2);
started = tic;
for bin = 1:180
    aa = a.Raw_Data(a.Angle>=bin-1 & a.Angle<bin);
    bb = b.Raw_Data(b.Angle>=bin-1 & b.Angle<bin);
    rng(1000+bin,'twister'); [p5000(bin),d] = permsim(aa,bb,5000);
    mc(bin,:) = d.MonteCarloCI95;
    if mod(bin,30)==0, fprintf('5,000-permutation audit: %d / 180 bins\n',bin); end
end
initialSeconds = toc(started);
[~,order] = sort(abs(p5000-0.05)); selected = sort(order(1:6));
p20000 = nan(180,1); pAuto = p20000; autoMC = nan(180,2); autoCount = p20000;
for bin = selected'
    aa = a.Raw_Data(a.Angle>=bin-1 & a.Angle<bin);
    bb = b.Raw_Data(b.Angle>=bin-1 & b.Angle<bin);
    rng(1000+bin,'twister'); p20000(bin) = permsim(aa,bb,20000);
    rng(1000+bin,'twister'); [pAuto(bin),d] = permsim(aa,bb,'auto', ...
        struct('TimeBudgetSeconds',300/180));
    autoCount(bin) = d.Iterations;
    autoMC(bin,:) = d.MonteCarloCI95;
    fprintf('Angle %.1f: 5k %.5f; 20k %.5f; auto %.5f [%.5f, %.5f]\n', ...
        angle(bin),p5000(bin),p20000(bin),pAuto(bin),autoMC(bin,:));
end
bin = 43;
aa = a.Raw_Data(a.Angle>=bin-1 & a.Angle<bin);
bb = b.Raw_Data(b.Angle>=bin-1 & b.Angle<bin);
oldTime = zeros(3,1); newTime = oldTime;
for repeat = 1:3
    rng(75); timer = tic; oldP = legacy(aa,bb,5000); oldTime(repeat) = toc(timer);
    rng(75); timer = tic; [newP,d] = permsim(aa,bb,5000); newTime(repeat) = toc(timer);
    assert(oldP==newP);
end
benchmark = struct('Pixels',numel(aa)+numel(bb),'Iterations',5000, ...
    'LegacySeconds',median(oldTime),'BatchedSeconds',median(newTime), ...
    'Speedup',median(oldTime)/median(newTime), ...
    'LegacyDataBytes',8*(numel(aa)+numel(bb))*5000, ...
    'BatchedDataBytes',8*(numel(aa)+numel(bb))*d.BatchSize);
report = table(angle,p5000,mc(:,1),mc(:,2),mc(:,1)<=0.05 & mc(:,2)>=0.05, ...
    p20000,pAuto,autoMC(:,1),autoMC(:,2),autoCount, ...
    'VariableNames',{'Angle','P5000','MC5000L','MC5000U','NearThreshold5000', ...
    'P20000','PAuto','MCAutoL','MCAutoU','AutoIterations'});
writetable(report,fullfile(folder,'permutation_precision_audit.csv'));
save(fullfile(folder,'permutation_precision_audit.mat'),'report','benchmark','plan','initialSeconds');
fig = figure('Visible','off','Color','w','Position',[100 100 1100 700]);
figureCleanup = onCleanup(@() close(fig));
tiles = tiledlayout(fig,2,1,'TileSpacing','compact');
ax = nexttile(tiles); errorbar(ax,angle,p5000,p5000-mc(:,1),mc(:,2)-p5000, ...
    '.','Color',[0.2 0.45 0.6],'CapSize',0); hold(ax,'on');
yline(ax,0.05,'--'); ylim(ax,[0 0.12]); xlim(ax,[0 180]); grid(ax,'on');
title(ax,'Current pixel test: 5,000 permutations, 1-degree bins');
xlabel(ax,'Angular distance (degrees)'); ylabel(ax,'Legacy p-value');
ax = nexttile(tiles); hold(ax,'on'); x = 1:numel(selected);
errorbar(ax,x-0.12,p5000(selected),p5000(selected)-mc(selected,1), ...
    mc(selected,2)-p5000(selected),'o','DisplayName','5,000');
errorbar(ax,x+0.12,pAuto(selected),pAuto(selected)-autoMC(selected,1), ...
    autoMC(selected,2)-pAuto(selected),'o','DisplayName','Auto: time-budgeted, 5k-20k');
yline(ax,0.05,'--','HandleVisibility','off'); grid(ax,'on');
xticks(ax,x); xticklabels(ax,string(angle(selected))); xlim(ax,[0.5 numel(selected)+0.5]);
xlabel(ax,'Angular position (degrees)'); ylabel(ax,'Legacy p-value');
title(ax,'Six bins closest to 0.05 in the initial run'); legend(ax,'Location','best');
title(tiles,'DMSO 4 dpi vs IWR1 4 dpi: permutation precision');
xlabel(tiles,'Bars: pointwise 95% Monte Carlo intervals, not biological confidence intervals.');
exportgraphics(fig,fullfile(folder,'permutation_precision_audit.png'),'Resolution',160);
disp(benchmark);
fprintf('5k run: %.1f s; %d / 180 intervals cross 0.05; %d zero estimates.\n', ...
    initialSeconds,nnz(report.NearThreshold5000),nnz(p5000==0));
end

function p = legacy(a,b,reps)
n = numel(a); m = numel(b); data = [a;b];
aa = zeros(n,reps); bb = zeros(m,reps);
for j = 1:reps
    order = randperm(n+m)'; aa(:,j) = data(order(1:n)); bb(:,j) = data(order(n+1:end));
end
p = nnz(abs(median(aa,1)-median(bb,1))>=abs(median(a)-median(b)))/reps;
end
