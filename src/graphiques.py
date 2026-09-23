"""
Relit resultats.hdf5 (produit par poisson_hdf5.py) et trace :
- la courbe de convergence de l'erreur L2 en fonction du pas h ;
- la solution approchée pour N = 50 en surface 3D.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

with h5py.File("../results/resultats.hdf5", "r") as f:
    liste_N = sorted(f[g].attrs["N"] for g in f.keys() if g.startswith("Sol_N"))
    erreurs = f["L2_error"][()]
    g50 = f["Sol_N50"]
    N50 = int(g50.attrs["N"])
    U50 = g50["Solution"][()]

# Convergence
h = np.array([1.0 / (N + 1) for N in liste_N])
pente = np.polyfit(np.log(h), np.log(erreurs), 1)[0]
plt.figure(figsize=(6, 4.5))
plt.loglog(h, erreurs, "o-", label=f"Erreur L2 (pente mesurée ≈ {pente:.2f})")
plt.loglog(h, erreurs[0] * (h / h[0]) ** 2, "--", color="gray", label="Pente 2 (référence)")
plt.xlabel("Pas d'espace h")
plt.ylabel("Erreur L2")
plt.title("Poisson 2D — convergence du schéma 5 points")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("../results/convergence.png", dpi=150)
print(f"Pente de convergence mesurée : {pente:.2f}")

# Surface de la solution (bords à 0 ajoutés)
M = N50 + 2
U = np.zeros((M, M))
U[1:-1, 1:-1] = U50.reshape((N50, N50))
x = np.linspace(0, 1, M)
X, Y = np.meshgrid(x, x)
fig = plt.figure(figsize=(6.5, 5))
ax = fig.add_subplot(projection="3d")
ax.plot_surface(X, Y, U, cmap="coolwarm", linewidth=0)
ax.set_title(f"Solution approchée, N = {N50}")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("u")
plt.tight_layout()
plt.savefig("../results/solution_N50.png", dpi=150)
print("Graphiques enregistrés dans results/")
