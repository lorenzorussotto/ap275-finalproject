import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd

# Read the CSV file
df = pd.read_csv("BaBiO3_ecut_optimization.csv")

# Plot nk on x-axis and energy on y-axis
plt.plot(df['nk'], df['energy'], marker='o', linestyle='-')
plt.xlabel('Plane wave energy cutoff (Ry)')
plt.ylabel('Energy (eV/atom)')
plt.title('Energy convergence vs. Plane wabe energy cutoff')
plt.grid(True)
plt.savefig('convergence_ecut.jpg')
plt.show()