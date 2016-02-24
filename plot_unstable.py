import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy import io

mpl.rcParams['lines.linewidth'] = 1.5

variables = io.loadmat('./results/2016-02-18-steady-piston-bc/unstable.mat')

grid = variables['grid'][0]
znd_all = variables['znd_all']

znd_u = znd_all['u'][0][0]
znd_l = znd_all['l'][0][0]
znd_w = znd_all['w'][0][0]

pert_real_u = variables['sol'][0, :]
pert_imag_u = variables['sol'][1, :]
pert_real_l = variables['sol'][2, :]
pert_imag_l = variables['sol'][3, :]

fig = plt.figure(figsize=(4.5, 8.3436))

ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(grid, znd_u, '-', label='$u$')
ax1.plot(grid, znd_l, '--', label='$\lambda$')
ax1.plot(grid, znd_w, '-.', label='$\omega$')
ax1.set_ylabel('ZND solution')
ax1.legend(loc='upper left')

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(grid, pert_real_u, label='Re $u$')
ax2.plot(grid, pert_real_l, label='Re $\lambda$')
ax2.set_ylabel('Real part')
ax2.legend(loc='lower left')

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(grid, pert_imag_u, label='Im $u$')
ax3.plot(grid, pert_imag_l, label='Im $\lambda$')
ax3.set_xlabel('$x$')
ax3.set_ylabel('Imaginary part')
ax3.legend(loc='upper left')

plt.tight_layout()
#plt.show()
plt.savefig('./results/2016-02-18-steady-piston-bc/unstable.pdf')
