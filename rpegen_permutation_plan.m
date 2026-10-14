function plan = rpegen_permutation_plan(halfWidth,referenceP)
% Fixed-run planning for Monte Carlo precision, not biological power.
% Normal-approximation 95% half-width at a chosen tail probability.
% Use referenceP=0.5 for worst-case absolute precision over all p-values.
if nargin < 1, halfWidth = 0.002; end
if nargin < 2, referenceP = 0.05; end
validateattributes(halfWidth,{'numeric'},{'scalar','real','finite','>',0,'<',1});
validateattributes(referenceP,{'numeric'},{'scalar','real','finite','>',0,'<',1});
reps = ceil(1.959963984540054^2*referenceP*(1-referenceP)/halfWidth^2);
plan = struct('Iterations',reps,'ReferenceP',referenceP, ...
    'ApproximateHalfWidth95',halfWidth,'Resolution',1/reps);
end
