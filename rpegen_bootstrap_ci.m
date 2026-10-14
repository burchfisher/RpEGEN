function result = rpegen_bootstrap_ci(raw,eyeNames,width,iterations,seed)
% Pointwise percentile intervals for pooled medians, resampling whole eyes.
% The same eye draws are used at every angle; pixels never cross angle bins.
if nargin < 4, iterations = 5000; end
if nargin < 5, seed = 20261012; end
validateattributes(width,{'numeric'},{'scalar','integer','positive'});
assert(ismember(width,[1 5 10]),'RpEGEN:Bootstrap','Unsupported angle bin width.');
validateattributes(iterations,{'numeric'},{'scalar','integer','>=',100});
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative','<=',2^32-1});
assert(istable(raw) && all(ismember({'IMG_name','Angle','Raw_Data'}, ...
    raw.Properties.VariableNames)),'RpEGEN:Bootstrap','Eye IDs and raw pixels are required.');
eyeNames = string(eyeNames(:)); n = numel(eyeNames);
assert(n >= 2 && ~any(ismissing(eyeNames) | strlength(eyeNames)==0) && ...
    numel(unique(eyeNames))==n,'RpEGEN:Bootstrap','At least two distinct eye IDs are required.');
[known,ids] = ismember(string(raw.IMG_name),eyeNames);
assert(all(known) && numel(unique(ids))==n,'RpEGEN:Bootstrap', ...
    'Raw pixel eye IDs must match the complete eye list.');
values = double(raw.Raw_Data); angles = double(raw.Angle);
assert(iscolumn(values) && iscolumn(angles) && all(isfinite(values)) && ...
    all(values>=0 & values<=255 & values==fix(values)) && ...
    all(isfinite(angles) & angles>=0 & angles<=180), ...
    'RpEGEN:Bootstrap','Expected finite 8-bit intensities and angles from 0 to 180.');
% A local stream makes display changes reproducible without changing permsim's RNG.
stream = RandStream('mt19937ar','Seed',seed);
draws = randi(stream,n,n,iterations);
drawNumber = repmat(1:iterations,n,1);
weights = accumarray([draws(:),drawNumber(:)],1,[n iterations]);
count = 180/width; lower = nan(count,1); upper = lower; pooled = lower;
contributing = zeros(count,1); emptyDraws = zeros(count,1);
for b = 1:count
    inBin = angles >= (b-1)*width & angles < b*width;
    histogram = accumarray([values(inBin)+1,ids(inBin)],1,[256 n]);
    contributing(b) = nnz(sum(histogram,1));
    pooled(b) = countMedian(sum(histogram,2));
    estimates = nan(1,iterations);
    for first = 1:512:iterations
        last = min(iterations,first+511);
        estimates(first:last) = countMedian(histogram*weights(:,first:last));
    end
    emptyDraws(b) = nnz(~isfinite(estimates));
    % Do not silently condition on nonempty draws or imply precision from one eye.
    if contributing(b) >= 2 && emptyDraws(b)==0
        limits = prctile(estimates,[2.5 97.5]);
        lower(b) = limits(1); upper(b) = limits(2);
    end
end
result = struct('Angle',((width/2):width:180)', 'Median',pooled, ...
    'Lower',lower,'Upper',upper,'ContributingEyes',contributing, ...
    'EmptyDraws',emptyDraws,'Iterations',iterations,'Seed',seed);
end

function value = countMedian(counts)
cumulative = cumsum(counts,1); total = cumulative(end,:);
value = (sum(cumulative < ceil(total/2),1)+ ...
    sum(cumulative < floor(total/2)+1,1))/2;
value(total==0) = NaN;
end
