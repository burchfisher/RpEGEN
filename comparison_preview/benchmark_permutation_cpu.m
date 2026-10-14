function benchmark_permutation_cpu
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
a = load(fullfile(root,'JoVE_dataset','Output','DMSO_4dpi.mat'));
b = load(fullfile(root,'JoVE_dataset','Output','IWR1_4dpi.mat'));
rawA = a.DMSO_4dpi.RAW_data; rawB = b.IWR1_4dpi.RAW_data;
before = rng;
[serial,serialInfo] = rpegen_permutation_profile(rawA,rawB,10,1000);
assert(isequal(before,rng));
[parallel,parallelInfo] = rpegen_permutation_profile(rawA,rawB,10,1000,struct('UseParallel',true));
assert(isequal(serial,parallel) && isequal(before,rng));
[parallel,parallelWarm] = rpegen_permutation_profile(rawA,rawB,10,1000,struct('UseParallel',true));
assert(isequal(serial,parallel));
speedup = serialInfo.ElapsedSeconds/parallelWarm.ElapsedSeconds;
fprintf('CPU benchmark: serial %.2f s, parallel cold %.2f s, warm %.2f s, speedup %.2fx\n', ...
    serialInfo.ElapsedSeconds,parallelInfo.ElapsedSeconds,parallelWarm.ElapsedSeconds,speedup);
save(fullfile(root,'comparison_preview','permutation_cpu_benchmark.mat'), ...
    'serialInfo','parallelInfo','parallelWarm','speedup');
% Confirm that auto mode works on thread workers and reports its actual count.
[automatic,autoInfo] = rpegen_permutation_profile(rawA,rawB,10,'auto', ...
    struct('UseParallel',true,'MinIterations',100,'MaxIterations',200,'TimeBudgetSeconds',30));
assert(all(automatic.Iterations>=100 & automatic.Iterations<=200));
assert(autoInfo.Parallel);
disp('Serial/parallel results, isolated RNG, and threaded auto mode passed.');
end
