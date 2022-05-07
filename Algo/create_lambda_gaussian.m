function lambda = create_lambda_gaussian(n, p, fdr)
  function o = w(k)
      o = 1 / max(1, n - k - 1);
  end
  lambda_bhq = create_lambda_bhq(n, p, fdr);
  lambda = zeros(p,1);
  lambda(1) = lambda_bhq(1);
  if (p >= 2)
    sum_sq = 0;
    for i= 2:p
      sum_sq = sum_sq + lambda(i-1)^2;
      lambda(i) = lambda_bhq(i) * sqrt(1 + w(i-1) * sum_sq);
    end
  end
end