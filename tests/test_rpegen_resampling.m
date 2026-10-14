function test_rpegen_resampling
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
previous = rng; cleanup = onCleanup(@() rng(previous));
for reps = [1 257 1025]
    for data = { {[1;2;2;7],[3;4;9]}, {ones(4,1),ones(6,1)}, ...
            {single([1;5]),single([3;7;10])}, {uint8([1;5]),uint8([3;7;10])} }
        pair = data{1}; rng(123); expected = legacy(pair{1},pair{2},reps); after = rng;
        rng(123); [actual,d] = permsim(pair{1},pair{2},reps);
        assert(isequal(actual,expected) && isequal(rng,after));
        assert(d.Exceedances/reps==actual);
        assert(d.MonteCarloCI95(1)<=actual && d.MonteCarloCI95(2)>=actual);
    end
end
plan = rpegen_permutation_plan; assert(plan.Iterations==45618);
assert(rpegen_permutation_plan(0.001).Iterations==182470);
rng(10); [p,d] = permsim([0;1],[2;3],'auto'); after = rng;
rng(10); q = permsim([0;1],[2;3],d.Iterations);
assert(p==q && d.Automatic && isequal(after,rng));
assert(d.Iterations>=10000 && d.Iterations<=50000);
assert(d.Planning.RequestedIterations==45618);
[~,d] = permsim([0;1],[2;3],'auto',struct('TimeBudgetSeconds',1e-10));
assert(d.Iterations==10000 && d.Planning.MinimumExceedsBudget);
[~,d] = permsim([0;1],[2;3],'auto',struct('HalfWidth',0.001,'TimeBudgetSeconds',60));
assert(d.Iterations==50000 && d.Planning.PrecisionTargetLimited);
[p,d] = permsim(ones(3,1),ones(5,1),100);
assert(p==1 && d.MonteCarloCI95(2)==1);
rng(3); [p,d] = permsim((0:99)',(100:199)',10);
assert(p==0 && d.ZeroExceedances && d.MonteCarloCI95(2)>0);
mustFail(@() permsim([],1,10)); mustFail(@() permsim([1;NaN],1,10));
mustFail(@() permsim([1;2],3,0)); mustFail(@() permsim([1;2],3,1.5));
mustFail(@() permsim([1;2],3,'auto',struct('HalfWidth',0)));
disp('Legacy seeded p-value and RNG parity, auto counts, and diagnostics passed.');

raw = table(["a";"a";"b";"b";"a";"a";"b";"b";"a"], ...
    [1;2;1;2;6;7;6;7;180],[0;2;10;12;100;102;200;202;255], ...
    'VariableNames',{'IMG_name','Angle','Raw_Data'});
before = rng; result = rpegen_bootstrap_ci(raw,["a";"b"],5,200,42);
assert(isequal(before,rng));
assert(isequaln(result,rpegen_bootstrap_ci(raw,["a";"b"],5,200,42)));
stream = RandStream('mt19937ar','Seed',42); draws = randi(stream,2,2,200);
for bin = 1:2
    values = nan(1,200);
    for draw = 1:200
        pixels = [];
        for eye = draws(:,draw)'
            names = ["a";"b"];
            rows = raw.IMG_name==names(eye) & raw.Angle>=(bin-1)*5 & raw.Angle<bin*5;
            pixels = [pixels;raw.Raw_Data(rows)]; %#ok<AGROW>
        end
        values(draw) = median(pixels);
    end
    assert(isequal(prctile(values,[2.5 97.5]),[result.Lower(bin),result.Upper(bin)]));
end
assert(all(isnan(result.Lower(3:end))) && isnan(result.Median(end)));
assert(result.Median(1)==6 && result.Median(2)==151);
sparse = raw; sparse(sparse.Angle<5 & sparse.IMG_name=="b",:) = [];
result = rpegen_bootstrap_ci(sparse,["a";"b"],5,200,42);
assert(isnan(result.Lower(1)) && result.EmptyDraws(1)>0);
mustFail(@() rpegen_bootstrap_ci(raw,["a";"a"],5));
disp('Whole-eye bootstrap matches literal pixel resampling; RNG and sparse bins passed.');
end

function mustFail(callback)
failed = false;
try, callback(); catch, failed = true; end
assert(failed,'Invalid input unexpectedly accepted.');
end

function p = legacy(a,b,reps)
n = numel(a); m = numel(b); data = [a;b];
aa = zeros(n,reps); bb = zeros(m,reps);
for j = 1:reps
    order = randperm(n+m)'; aa(:,j) = data(order(1:n)); bb(:,j) = data(order(n+1:end));
end
p = nnz(abs(median(aa,1)-median(bb,1))>=abs(median(a)-median(b)))/reps;
end
