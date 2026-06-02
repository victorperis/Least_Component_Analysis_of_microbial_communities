import numpy as np
import matplotlib.pyplot as plt


def deltaL(lam):
    out = 0.5*(lam-1-np.log(lam))
    return out

lam = np.linspace(1e-10,6,1000)

fig,ax = plt.subplots(figsize=(4,2))
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
ax.plot(lam,deltaL(lam),
        'k-',
        linewidth=3)
ax.set_xlabel(r'$\lambda_\alpha$')
ax.set_ylabel(r'$\Delta \mathcal{L} (\lambda_\alpha)$')
ax.set_ylim((-0.1,3))
ax.set_xticks([0,1,2,4,6])
ax.set_xticklabels([0,1,2,4,6])
plt.tight_layout()
plt.savefig('pattern_loglikelihood.png',format='png',dpi=300,transparent=True)
plt.show()
