% #' Lambda sequences for SLOPE
% #'
% #' Computes \eqn{\lambda} sequences for SLOPE according to several pre-defined methods.
% #'
% #' @param n number of observations
% #' @param p number of variables
% #' @param fdr target False Discovery Rate (FDR)
% #' @param method method to use for computing \eqn{\lambda} (see Details)
% #'
% #' @details The following methods for computing \eqn{\lambda} are supported:
% #' \itemize{
% #'  \item \code{bhq}: Computes sequence inspired by Benjamini-Hochberg (BHq)
% #'    procedure
% #'  \item \code{gaussian}: Computes modified BHq sequence inspired by
% #'    Gaussian designs
% #' }
% #'
% #' @rdname lambda
% #' @export
function out = create_lambda(n, p, fdr, method)

  switch (method)
      case 'bhq'
          out = create_lambda_bhq(n, p, fdr);
      case 'gaussian'
          out = create_lambda_gaussian_truncated(n, p, fdr);
  end
end

function out = create_lambda_bhq(n, p, fdr) 
  q = [1:p] * fdr / (2*p);
  out = norminv(1 - q);
end

function lambda = create_lambda_gaussian(n, p, fdr)
  function o = w(k)
      o = 1 / max(1, n - k - 1);
  end
  lambda_bhq = create_lambda_bhq(n, p, fdr);
  lambda = rep(0,p);
  lambda(1) <- lambda_bhq(1);
  if (p >= 2)
    sum_sq = 0;
    for i= 2:p
      sum_sq = sum_sq + lambda(i-1)^2;
      lambda(i) = lambda_bhq(i) * sqrt(1 + w(i-1) * sum_sq);
    end
  end
end

function lambda = create_lambda_gaussian_truncated(n, p, fdr)
  lambda = create_lambda_gaussian(n, p, fdr);
  [minlam, k] = min(lambda);
  lambda(k:p) = lambda(k);
end