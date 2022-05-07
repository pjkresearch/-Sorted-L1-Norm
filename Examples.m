%% Index Tracking
clear; close all; clc;
addpath('Algo');


%Set Parameters
num_lambdas  = 100; % number of lambdas
q            = 0.1; % determines how fast lambda decreases in SL1
match_idx    = 1; % match lambda_slope(match_idx) = lambda_lasso
phi          = 2; %as in problem definition 

options = struct('max_iter', 8000, ...
                 'tolInfeas', 1e-7, ...
                 'tolRelGap', 1e-7); 
             

% Sample asset return data
% Create Matrices
k = 3;
n = 50;
p = 12;
NumberBlocks = 3;
f1=round(p/NumberBlocks);  

C = [ repmat([0.77 0.64     0], f1, 1);...
   repmat([0.9  0    -0.42], f1, 1);...
   repmat([0    0.31  0.64], f1,1)]';

B = mvnrnd(zeros(1,k), eye(3), n);
SigmaEpsilon = 0.05;
Epsilon = mvnrnd(zeros(1,p), repmat(SigmaEpsilon, 1, p), n);
X = B*C + Epsilon;

beta = 0.1*ones(p,1);
Y = X*beta + randn(n,1);

[T, n] = size(X);


VCM_true = (C' * cov(B))* C + diag(diag(cov(Epsilon))); 
ISexpCov = cov(X); 
ISexpRet = mean(X);

%Set equally weighted portfolio as starting point
e        = ones(1,n);
eq_port  = e'*(1/n);



%% Index Tracking: Slope and Budget Constraint
disp('Compute SLOPE for Index Tracking')
for i =1:num_lambdas
    lambda_min   =    -.5;
    lambda_max   =      1;   
    lambda_lasso = logspace(lambda_min, lambda_max, num_lambdas);

    beta_hat=[];
    lambda_slope = create_lambda(T, n, q, 'bhq');
    lambda_slope = lambda_lasso(i).* lambda_slope;

    if i > 1
    soln = regADM_bd(X, Y, 'SL1', lambda_slope, zeros(n,1), ones(n,1), false, options, soln.w);
    else
    soln = regADM_bd(X, Y, 'SL1', lambda_slope, zeros(n,1), ones(n,1), false, options, eq_port); 
    end
    
    beta_hat = soln.w;
    t=tabulate(abs(beta_hat));
    fprintf('Sample - i=%d lam=%5.3e nnz=%d npar=%d n_short=%d status=%s\n', i, lambda_lasso(i), nnz(beta_hat), size(t,1), nnz(beta_hat<0), soln.status);

    
    %Compute Statistics
    index_clumps_temp(i)                = length(t(:,2));
    index_clumps(i)                     = max([nnz((t(:,2)>1)) 1]);
    index_beta_objective_lasso(i)       = soln.obj; 
    index_nactive(i)              = nnz(beta_hat);
    index_short(i)                = nnz(beta_hat<0);
    index_RSS(i)                  = sum((Y - X*beta_hat).^2);
    index_BIC(i)                  = (-2*log(index_RSS(i)/T))+ log(T)*(index_nactive(i));      % Bayesian information criterion       
    index_tot_betahat(:,i)   = beta_hat;

end
   

%% Mean-Variance Optimization

for i =1:num_lambdas
    disp('Compute SLOPE for Mean-Variance Optimization')
    lambda_min   =     -2;
    lambda_max   =      1;   
    lambda_lasso = logspace(lambda_min, lambda_max, num_lambdas);

    beta_hat=[];
    lambda_slope = create_lambda(T, n, q, 'bhq');
    lambda_slope = (lambda_lasso(i) / lambda_slope(match_idx)) .* lambda_slope;

    % Compute Optimal Weights using Estimated CoVar
    if(i>1)
        soln=varADM(ISexpCov, ISexpRet, phi, 'SL1', lambda_slope, false, options, soln.w);
    else
        soln=varADM(ISexpCov, ISexpRet, phi, 'SL1', lambda_slope,false, options);
    end
    
    beta_hat = soln.w;
    t=tabulate(abs(beta_hat));
    fprintf('Sample - i=%d lam=%5.3e nnz=%d npar=%d n_short=%d status=%s\n', i, lambda_lasso(i), nnz(beta_hat), size(t,1), nnz(beta_hat<0), soln.status);

    beta_hat = soln.w;
    
    t1=tabulate(abs(beta_hat));
    
    %Compute Statistics
    port_weights(:, i)         = beta_hat;
    port_aw(i)                 = nnz(beta_hat);
    port_num_short_pos(i)        = nnz(beta_hat<0);
    port_sparsity(i)           = sum(abs(beta_hat)~=0);
    port_risk(i)               = sqrt((beta_hat'*ISexpCov)*beta_hat);
end








