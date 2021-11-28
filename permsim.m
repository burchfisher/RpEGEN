function p = permsim(dataA,dataB,reps)
% Created by G. Burch Fisher beginning on 11/21/21
% Last updated 11/24/21

% Calculates a p-value based on comparing the medians of a permutation simulation between two datasets
% Can be easily modified to use mean or any other test statistic by editing below.

% Concatenate the two arrays dataA and dataB into data
lenA = length(dataA);
lenB = length(dataB);
data = [dataA; dataB];

% Generate permutations equal to the number of repetitions
permA_data = zeros(lenA,reps);
permB_data = zeros(lenB,reps);

for x = 1:reps;
    perm = transpose(randperm(lenA + lenB));
    permA_data(:,x) = data(perm(1:lenA));
    permB_data(:,x) = data(perm(lenA+1:end));
end

% Calculate the difference in medians for each of the datasets
samples = median(permA_data,1) - median(permB_data,1);

% Calculate the test statistic and p-value
test_stat = median(dataA)-median(dataB);
p = sum(abs(samples) >= abs(test_stat))/reps;
end