import numpy as np
import scipy.stats as stats

# observed contingency table in the slides
obs1 = np.array([[28, 42, 56],[95, 31, 17]])
obs2 = np.array([[448,482,466],[435,511,497]])

obs = obs1 + obs2
# Pearson's chi-squared test (two-tailed)
print('>>> Pearson\'s chi-squared test <<<')
chi, p, dof, exp = stats.chi2_contingency(obs)
print(' chi-squared statistic: %.2f' % chi)
print(' degree of freedom: %d' % dof)
print(' p-value: %.2e' % p)



# Armitage test for linear trend (two-tailed)
print('\n>>> Armitage test for linear trend <<<')
C, R = np.sum(obs, axis=0), np.sum(obs, axis=1)
N = np.sum(C)
t = np.array([0, 1, 2])	# linear weights
T = np.sum(t*(obs[0]*R[1]-obs[1]*R[0]))
var =  R[0]*R[1]/N * (np.sum(t**2*C*(N-C)) - 2*t[1]*t[2]*C[1]*C[2])
z = T/var**0.5
p = 2*stats.norm.sf(abs(z))
print(' z-statistic: %.2f' % z)
print(' p-value: %.2e' % p)

