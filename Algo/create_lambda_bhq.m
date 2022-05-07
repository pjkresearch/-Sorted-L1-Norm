function out = create_lambda_bhq(n, p, fdr) 
  q = [1:p] * fdr / (2*p);
  out = norminv(1 - q);
end