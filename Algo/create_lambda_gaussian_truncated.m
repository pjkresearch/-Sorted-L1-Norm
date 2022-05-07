function lambda = create_lambda_gaussian_truncated(n, p, fdr)
  lambda = create_lambda_gaussian(n, p, fdr);
  [minlam, k] = min(lambda);
  lambda(k:p) = lambda(k);
end