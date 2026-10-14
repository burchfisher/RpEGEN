function [p,diagnostics] = permsim(dataA,dataB,reps,options)
% Created by G. Burch Fisher beginning on 11/21/21
% Last updated 11/24/21

% Calculates a p-value based on comparing the medians of a permutation simulation between two datasets
% Can be easily modified to use mean or any other test statistic by editing below.

% The original pixel shuffle, two-sided statistic, and b/reps formula are
% unchanged. The optional second output describes Monte Carlo precision only.
% reps='auto' benchmarks first, then fixes a count without looking at p.
% options: HalfWidth (0.002), TimeBudgetSeconds (2), MinIterations (10000),
% MaxIterations (50000). The runtime budget is approximate; the minimum wins.
% Numeric reps keeps manual counts and the original RNG sequence exactly.
if nargin < 3, reps = 'auto'; end
if nargin < 4, options = struct; end
cancelFcn = [];
if isstruct(options) && isfield(options,'CancelFcn')
    cancelFcn = options.CancelFcn;
    assert(isempty(cancelFcn) || isa(cancelFcn,'function_handle'), ...
        'RpEGEN:Permutation','CancelFcn must be a function handle.');
end
validateattributes(dataA,{'numeric'},{'column','nonempty','real','finite'});
validateattributes(dataB,{'numeric'},{'column','nonempty','real','finite'});
automatic = (ischar(reps) || isstring(reps)) && strcmp(reps,'auto');
if automatic
    settings = inputParser;
    addParameter(settings,'HalfWidth',0.002,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<1);
    addParameter(settings,'TimeBudgetSeconds',2,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
    addParameter(settings,'MinIterations',10000,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=100 && x==fix(x));
    addParameter(settings,'MaxIterations',50000,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=100 && x==fix(x));
    addParameter(settings,'CancelFcn',[],@(x) isempty(x) || isa(x,'function_handle'));
    parse(settings,options); config = settings.Results;
    assert(config.MaxIterations>=config.MinIterations,'RpEGEN:Permutation', ...
        'MaxIterations must be at least MinIterations.');
    plan = rpegen_permutation_plan(config.HalfWidth,0.05);
    % Discard calibration results and restore RNG: count selection never sees p.
    previous = rng; restore = onCleanup(@() rng(previous));
    permsim(dataA,dataB,32);
    calibration = tic; permsim(dataA,dataB,128); secondsPerDraw = toc(calibration)/128;
    clear restore;
    desired = min(plan.Iterations,config.MaxIterations);
    affordable = floor(config.TimeBudgetSeconds/max(secondsPerDraw,eps));
    reps = max(config.MinIterations,min(desired,affordable));
    plan.RequestedIterations = plan.Iterations; plan.Iterations = reps;
    plan.TimeBudgetSeconds = config.TimeBudgetSeconds;
    plan.EstimatedSeconds = reps*secondsPerDraw;
    plan.PrecisionTargetLimited = reps < plan.RequestedIterations;
    plan.MinimumExceedsBudget = affordable < config.MinIterations;
    plan.ApproximateHalfWidth95 = 1.959963984540054*sqrt(0.05*0.95/reps);
    plan.Resolution = 1/reps;
end
validateattributes(reps,{'numeric'},{'scalar','integer','positive','finite'});
lenA = length(dataA);
lenB = length(dataB);
data = [dataA; dataB];
test_stat = median(dataA)-median(dataB);
% Bound the two shuffled data arrays to about 32 MiB, independent of reps.
batchSize = max(1,min(256,floor(32*1024^2/(8*(lenA+lenB)))));
exceedances = 0;
for first = 1:batchSize:reps
    if ~isempty(cancelFcn) && cancelFcn()
        error('RpEGEN:Cancelled','Permutation test cancelled.');
    end
    batch = min(batchSize,reps-first+1);
    permA_data = zeros(lenA,batch);
    permB_data = zeros(lenB,batch);
    for x = 1:batch
        perm = transpose(randperm(lenA + lenB));
        permA_data(:,x) = data(perm(1:lenA));
        permB_data(:,x) = data(perm(lenA+1:end));
    end
    samples = median(permA_data,1) - median(permB_data,1);
    exceedances = exceedances + nnz(abs(samples) >= abs(test_stat));
end
p = exceedances/reps;
if nargout > 1
    [~,limits] = binofit(exceedances,reps,0.05);
    diagnostics = struct('Iterations',reps,'Exceedances',exceedances, ...
        'MonteCarloCI95',limits,'Resolution',1/reps, ...
        'ZeroExceedances',exceedances==0, ...
        'NearThreshold',limits(1)<=0.05 && limits(2)>=0.05, ...
        'BatchSize',batchSize,'TestStatistic',test_stat,'Automatic',automatic);
    if automatic, diagnostics.Planning = plan; end
end
end
