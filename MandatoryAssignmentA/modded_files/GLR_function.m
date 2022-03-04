clear all;
clc;

load('classroomfile_7_2.mat');
M = 25;

mu_0 = 0;
sigma = 0.1;

g = zeros(size(r));     % Test statistics
idx = ones(size(r));    % Index of fault occurence sample estimation

%sum = zeros(M,1):

for k = M:length(r)

    g(k) = 1/(2*sigma^2*M)*max(sum(r(k-M+1:k)-mu_0))^2;
    sum(r(k-M+1:k)-mu_0)^2;
    
end

