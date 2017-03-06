lt = function(x) x[lower.tri(x)]
ztrans = function(x) 0.5 * log((1 + x) / (1 - x))
inv_ztrans = function(x){ (exp(2*x) - 1) / (exp(2*x) + 1) }
