function out = create_lambda(n, p, fdr, method)

  switch (method)
      case 'bhq'
          out = create_lambda_bhq(n, p, fdr);
      case 'gaussian'
          out = create_lambda_gaussian_truncated(n, p, fdr);
  end
end