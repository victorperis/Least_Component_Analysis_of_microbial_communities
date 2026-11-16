from math import comb
from fractions import Fraction
from scipy.stats import hypergeom  # optional cross-check

# Values from Methods 4.4
M, K, m, kstar = 276, 73, 99, 72  # samples, taxon A count, taxon B count, observed co-occurrences

def P(k):
    """Exact hypergeometric pmf, P(k) = C(K,k) C(M-K,m-k) / C(M,m)."""
    return Fraction(comb(K, k) * comb(M - K, m - k), comb(M, m))

k_max = min(K, m)   # largest possible overlap (= 73 here)

# Sanity check: pmf sums to 1 (exact arithmetic)
assert sum(P(k) for k in range(max(0, K + m - M), k_max + 1)) == 1

p_ge = sum(P(k) for k in range(kstar,     k_max + 1))   # sum_{k >= k*} (what I think was intended)

print(f"P(k >= {kstar}) = {float(p_ge):.3e}")

# Cross-check with scipy (arguments: population size, successes in population, draws)
rv = hypergeom(M, K, m)
print(f"scipy sf({kstar-1}) = {rv.sf(kstar-1):.3e}   (this is P(k >= k*))")
