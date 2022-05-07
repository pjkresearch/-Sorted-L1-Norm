function soln = regADM_bd(X, Y, pen, lambda, lb, ub, verbose, options, init_w)
% regADM     A solver for penalized regression with an equality constraint
%
% soln = regADM(X, Y, lambda, options)
%
% The solver finds the solution of 
%
%    min_{w in R^K}   .5 || Y - Xw ||^2 + pen_lambda(w)
%               s.t.  sum_i w_i = 1
%
% Other arguments are:
%
%    pen          The type of penality term. 'L1' or 'SL1'
%
%    lambda       A scaler (if pen='L1') or a vector of length K 
%                 (if pen='SL1') sorted in a decreasing order: 
%                     lambdas_1 >= lambdas_2 >= ... >= lambdas_K >= 0
%
%    lb           A vector of lower bounds of the weights
%
%    ub           A vector of upper bounds of the weights
%
%    verbose      [Optional: default=false] Verbose output of algorithm
%                 progress
%
%    options      [Optional] a structure to control optimization
%                 .max_iter     Max no of iterations [default: 10000]
%                 .tolInfeas    Infeasibility tolerance [default: 1e-6]
%                 .tolRelGap    Relative primal-dual gap [default: 1e-6]
%                 .rho          A parameter for ADM [default: normest(X)]
%
%    init_w       [Optional] initial point of w. Default is a vector with
%                 all element = 1/K.
%
% The output soln is a structure with fields:
%   .obj          Optimal objective value
%   .w            Optimal solution
%   .status       Solver status

% Author: Sangkyun Lee <sangkyun.lee@cs.tu-dortmund.de>
% Date: Oct 7, 2016
%

    % Constants for exit status
    STATUS_RUNNING    = 0;
    STATUS_OPTIMAL    = 1;
    STATUS_ITERATIONS = 2;
    STATUS_MSG = {'Optimal','Iteration limit reached'};

    K = size(X,2);
    if(size(Y,2) ~= 1)
        Y = Y';
    end
    
    switch pen
        case 'L1'
            if(verbose)
                fprintf('L1 penalization\n\n');
            end
            if(length(lambda) > 1)
                fprintf('> Error! lambda should be a scalar\n\n');
                return;
            end
        case 'SL1'
            if(verbose)
                fprintf('SL1 penalization\n\n');
            end
            if(length(lambda) ~= K)
                fprintf('> Error! lambda should be a vector of length %d, now %d\n\n', K, length(lambda));
                return;
            end
            if(size(lambda,2) ~= 1)
                lambda = lambda';
            end
        otherwise
            fprintf('> Error! pen must be either ''L1'' or ''SL1''\n\n');
            return;
    end
    
    % optional arguments
    if ~exist('verbose','var') || isempty(verbose)
        verbose = false;
    end
    if ~exist('options','var') || isempty(options)
        options = struct(); 
    end
    if ~exist('init_w','var') || isempty(init_w)
        init_w = ones(K,1)/K;
    end
    if length(lb) ~= K || length(ub) ~= K
        fprintf('> Error! lb and ub should have the length of weights %d\n\n', K)
        return;
    end

    max_iter = getDefaultField(options,'max_iter',10000);
    tolInfeas  = getDefaultField(options,'tolInfeas',1e-6);
    tolRelGap  = getDefaultField(options,'tolRelGap',1e-6);
    rho0 = getDefaultField(options,'rho',normest(X));
    rho = rho0;

    % initiate variables
    w = init_w;
    v = w;
    alp = zeros(K,1);
    beta = 0;
    e = ones(K,1);
    
    status = STATUS_RUNNING;
    
    XX = X'*X;
    XY = X'*Y;
    
    [Q,R] = qr(XX + rho*(eye(K) + e*e'));
    
    if(verbose)
        fprintf('Iter:  Primal       Dual       PDGap        InfeasP     InfeasD     rho    \n');
        fprintf('---------------------------------------------------------------------------\n');
    end
    
    for j=1:max_iter
       
       w = clip(R\(Q'*(XY - alp - beta + rho*(v+1))), lb, ub);
       
       switch pen
           case 'L1'
           v = proxL1(w + (1/rho)*alp, lambda/rho);
           
           case 'SL1'
           v = proxSortedL1(w + (1/rho)*alp, lambda/rho);
       end
       
       alp = alp + rho*(w - v);
       beta = beta + rho*(sum(w) - 1); 

       Xw = X*w;
       r = Y - Xw;
       f = .5*(r'*r);
       switch pen
           case 'L1'
           obj_p = f + lambda*sum(abs(w));

           case 'SL1'
           sw = sort(abs(w), 'descend');
           obj_p = f  + lambda'*sw;
       end

       obj_d = f + (alp + beta)'*w - beta;
       infeas_p = max( norm(w-v,Inf), abs(sum(w) - 1) );

       switch pen
           case 'L1'
           infeas_d = max(max(abs(alp))-lambda, 0);

           case 'SL1'
           sa = sort(abs(alp), 'descend');
           infeas_d = max( max(max(cumsum(sa - lambda)), 0) );
       end
       pdgap = abs(obj_p - obj_d) / max(1, obj_p);

       if(verbose)
          fprintf('%04d: %5.3e   %5.3e   %5.3e    %5.3e   %5.3e   %5.3e\n', j, obj_p, obj_d, pdgap, infeas_p, infeas_d, rho);
       end

       if(pdgap < tolRelGap && max(infeas_p, infeas_d) < tolInfeas)
           if(verbose)
            fprintf('Optimal\n'); 
           end
           status = STATUS_OPTIMAL;
           break;
       end
       
    end
    
    if(status ~= STATUS_OPTIMAL)
        status = STATUS_ITERATIONS;
    end
    
    soln.w = v;
    soln.pen = pen;
    soln.alp = alp;
    soln.obj = obj_p;
    soln.status = STATUS_MSG{status};
    switch pen
       case 'L1'
       soln.J = lambda*sum(abs(w));

       case 'SL1'
       sw = sort(abs(v), 'descend');
       soln.J = lambda'*sw;
    end
    
end


% ------------------------------------------------------------------------
function opt = getDefaultField(data,field,default)
% ------------------------------------------------------------------------
   if isfield(data,field)
      opt = data.(field);
   else
      opt = default;
   end
end

% ------------------------------------------------------------------------
function x = proxL1(y,lambda)
% ------------------------------------------------------------------------
   % Normalization
   %y    = y(:);
   %sgn  = sign(y);
   x    = sign(y) .* max(abs(y) - lambda,0);
end

% ------------------------------------------------------------------------
function z = clip(x, lb, ub)
% ------------------------------------------------------------------------
   z    = max(lb, min(ub, x));
end