import pickle
import pylab as plt
import numpy as np

with open("tov_BSk19_RK4.pickle", "rb") as fr:
    data = pickle.load(fr)

data = np.matrix(data)
print(data)

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$R$ vs $\rho$", fontsize=16)
plt.xlabel(r'$r$', fontsize=14)
plt.ylabel(r'$\rho$', fontsize=14)

r = data[:,0]
m = data[:,1]
rho = data[:,2]

plt.plot(r, rho)

plt.legend(fontsize=12)
plt.grid()
plt.savefig("r_vs_rho.png")

# Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$R$ vs $m$", fontsize=16)
plt.xlabel(r'$r$', fontsize=14)
plt.ylabel(r'$m$', fontsize=14)

plt.plot(r, m)

plt.legend(fontsize=12)
plt.grid()
plt.savefig("r_vs_m.png")

p = 30000 * np.squeeze(np.asarray(rho)) ** 2

# Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$R$ vs $P$", fontsize=16)
plt.xlabel(r'$r$', fontsize=14)
plt.ylabel(r'$p$', fontsize=14)

plt.plot(r, p)

plt.legend(fontsize=12)
plt.grid()
plt.savefig("r_vs_p.png")
