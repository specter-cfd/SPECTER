import numpy as np
import matplotlib.pyplot as plt

# Plots variables in balance.txt as a function of time
# Assumes balance.txt is the output of an HD run
# To run execute in a shell 'python plot_balance.py'
# Execute in ipython with '%run plot_balance.py'

# Path to the data
path = '../../bin/'

# Reads balance.txt
#  t   = time
#  ene = energy (v^2)
#  ens = enstrophy (w^2)
#  eps = bulk energy injection rate
t, ene, ens, eps = np.loadtxt(path+'balance.txt', unpack=True)

# Plots energy vs. time in a window
# Change ene to ens or eps to plot enstrophy or
# injection rate instead, respectively.
plot.figure()
plt.plot(t, ene)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.show()

# Saves plot to an EPS file
plt.savefig('balance.eps', format='eps', dpi=600)
