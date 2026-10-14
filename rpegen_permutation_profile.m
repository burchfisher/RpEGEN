function [report,runInfo] = rpegen_permutation_profile(rawA,rawB,width,reps,options)
% Compare angle bins using the legacy pixel test, optionally on CPU workers.
% Random streams are fixed per bin, independent of worker scheduling.
if nargin < 3, width = 1; end
if nargin < 4, reps = 'auto'; end
if nargin < 5, options = struct; end
settings = inputParser;
addParameter(settings,'TimeBudgetSeconds',300,@positiveScalar);
addParameter(settings,'HalfWidth',0.002,@positiveScalar);
addParameter(settings,'MinIterations',10000,@positiveInteger);
addParameter(settings,'MaxIterations',50000,@positiveInteger);
addParameter(settings,'UseParallel',false,@(x) islogical(x) && isscalar(x));
addParameter(settings,'NumWorkers',4,@positiveInteger);
addParameter(settings,'ProgressFcn',[],@(x) isempty(x) || isa(x,'function_handle'));
addParameter(settings,'CancelFcn',[],@(x) isempty(x) || isa(x,'function_handle'));
addParameter(settings,'Seed',20261012,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0 && x==fix(x) && x<=2^32-181);
parse(settings,options); config = settings.Results;
automatic = (ischar(reps) || isstring(reps)) && isscalar(string(reps)) && strcmp(reps,'auto');
if ~automatic, validateattributes(reps,{'numeric'},{'scalar','integer','positive','finite'}); end
validateattributes(width,{'numeric'},{'scalar','integer','positive'});
assert(mod(180,width)==0,'RpEGEN:Permutation','Bin width must divide 180.');
assert(config.HalfWidth<1 && config.MinIterations>=100 && ...
    config.MaxIterations>=config.MinIterations,'RpEGEN:Permutation','Invalid precision or iteration limits.');
for raw = {rawA,rawB}
    entry = raw{1};
    assert(istable(entry) && all(ismember({'Angle','Raw_Data'},entry.Properties.VariableNames)), ...
        'RpEGEN:Permutation','Expected Angle and Raw_Data columns.');
    validateattributes(entry.Angle,{'numeric'},{'column','finite','real','>=',0,'<=',180});
    validateattributes(entry.Raw_Data,{'numeric'},{'column','finite','real'});
end
count = 180/width; a = cell(count,1); b = a;
progress(config,0,count,'Checking bins');
for bin = 1:count
    a{bin} = rawA.Raw_Data(rawA.Angle>=(bin-1)*width & rawA.Angle<bin*width);
    b{bin} = rawB.Raw_Data(rawB.Angle>=(bin-1)*width & rawB.Angle<bin*width);
    assert(~isempty(a{bin}) && ~isempty(b{bin}),'RpEGEN:Permutation', ...
        'Both groups need pixels in bin [%g, %g).',(bin-1)*width,bin*width);
end
started = tic; workers = 1; useParallel = false;
if config.UseParallel
    progress(config,0,count,'Starting CPU workers');
    if license('test','Distrib_Computing_Toolbox')
        try
            pool = gcp('nocreate');
            if isempty(pool), pool = parpool('Threads',config.NumWorkers); end
            workers = min(config.NumWorkers,pool.NumWorkers); useParallel = workers>1;
        catch problem
            warning('RpEGEN:ParallelFallback','Using serial execution: %s',problem.message);
        end
    else
        warning('RpEGEN:ParallelFallback','Parallel toolbox unavailable; using serial execution.');
    end
end
budget = max(eps,config.TimeBudgetSeconds-toc(started));
binOptions = struct('HalfWidth',config.HalfWidth,'MinIterations',config.MinIterations, ...
    'MaxIterations',config.MaxIterations,'TimeBudgetSeconds',budget*workers/count);
p = zeros(count,1); diagnostics = cell(count,1); seed = config.Seed;
interactive = ~isempty(config.ProgressFcn) || ~isempty(config.CancelFcn);
execution = 'serial';
if useParallel && interactive
    % Cancellable futures keep GUI callbacks on the client. The per-bin test
    % and seeds are identical to the noninteractive parfor path below.
    execution = 'parfeval';
    futures = parallel.FevalFuture.empty(0,1);
    cleanups = cell(count,1);
    for bin = 1:min(workers,count)
        future = parfeval(pool,@oneBin,2,a{bin},b{bin},reps,binOptions,seed+bin);
        futures(bin,1) = future;
        cleanups{bin} = onCleanup(@() cancel(future));
    end
    next = numel(futures)+1;
    completed = 0;
    while completed < count
        progress(config,completed,count,'Permuting pixels');
        [bin,value,detail] = fetchNext(futures,0.1);
        if ~isempty(bin)
            p(bin) = value; diagnostics{bin} = detail; completed = completed+1;
            if next <= count
                future = parfeval(pool,@oneBin,2,a{next},b{next},reps,binOptions,seed+next);
                futures(next,1) = future;
                cleanups{next} = onCleanup(@() cancel(future));
                next = next+1;
            end
        end
    end
elseif useParallel
    execution = 'parfor';
    parfor (bin = 1:count,workers)
        [p(bin),diagnostics{bin}] = oneBin(a{bin},b{bin},reps,binOptions,seed+bin);
    end
else
    binOptions.CancelFcn = config.CancelFcn;
    for bin = 1:count
        progress(config,bin-1,count,'Permuting pixels');
        [p(bin),diagnostics{bin}] = oneBin(a{bin},b{bin},reps,binOptions,seed+bin);
    end
end
progress(config,count,count,'Complete');
d = vertcat(diagnostics{:}); mc = vertcat(d.MonteCarloCI95);
report = table(((width/2):width:180)',p,[d.Iterations]',[d.Exceedances]', ...
    mc(:,1),mc(:,2),[d.NearThreshold]', ...
    'VariableNames',{'Angle','LegacyP','Iterations','Exceedances','MC95L','MC95U','NearThreshold'});
runInfo = struct('ElapsedSeconds',toc(started),'TimeBudgetSeconds',config.TimeBudgetSeconds, ...
    'Parallel',useParallel,'Workers',workers,'Seed',seed,'BinWidth',width,'Execution',execution);
end

function progress(config,done,total,message)
if ~isempty(config.ProgressFcn), config.ProgressFcn(done,total,message); end
if ~isempty(config.CancelFcn) && config.CancelFcn()
    error('RpEGEN:Cancelled','Permutation test cancelled.');
end
end

function [p,d] = oneBin(a,b,reps,options,seed)
previous = rng; cleanup = onCleanup(@() rng(previous));
rng(seed,'twister'); [p,d] = permsim(a,b,reps,options);
end

function valid = positiveScalar(value)
valid = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value) && value>0;
end

function valid = positiveInteger(value)
valid = positiveScalar(value) && value==fix(value);
end
