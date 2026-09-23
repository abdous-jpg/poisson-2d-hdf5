import numpy as np
import scipy.linalg as sc
import time
import h5py

#  Problème de Poisson 2D 
def f_source(x, y):
    return np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)

def u_analytique(x, y):
    return (1.0 / (8.0 * np.pi**2)) * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)

def solve_poisson(N):
    h = 1.0 / (N + 1)
    n2 = N * N
    
    # Grille interne
    pts = np.linspace(h, 1-h, N)
    X, Y = np.meshgrid(pts, pts)
    
    # Construction de la matrice (Laplacien 5 points)
    A = np.zeros((n2, n2))
    for j in range(N):
        for i in range(N):
            k = i + j * N
            A[k, k] = 4.0 / h**2
            # Voisins (i-1, i+1, j-1, j+1)
            for di, dj in [(-1,0), (1,0), (0,-1), (0,1)]:
                ni, nj = i + di, j + dj
                if 0 <= ni < N and 0 <= nj < N:
                    A[k, ni + nj * N] = -1.0 / h**2

    # Second membre
    b = f_source(X.ravel(), Y.ravel())

    # Résolution LU
    lu, piv = sc.lu_factor(A)
    u_approx = sc.lu_solve((lu, piv), b)
    
    return u_approx, h

#  Script Principal 
liste_N = [10, 20, 50, 80]
tab_err = []
tab_cpu = []

# Création du fichier HDF5
with h5py.File("../results/resultats.hdf5", "w") as fic:
    for N in liste_N:
        t_start = time.time()
        U, h = solve_poisson(N)
        t_end = time.time() - t_start
        
        # Calcul de l'erreur L2
        x_c = np.linspace(h, 1-h, N)
        xc, yc = np.meshgrid(x_c, x_c)
        uex = u_analytique(xc, yc).ravel()
        erreur = h * np.linalg.norm(U - uex)
        
        tab_err.append(erreur)
        tab_cpu.append(t_end)
        
        # Organisation en groupes selon le TD
        nom_groupe = "Sol_N" + str(N)
        gp = fic.create_group(nom_groupe)
        gp.create_dataset("Solution", data=U)
        gp.attrs['N'] = N 
        
        print(f"N = {N} calculé. Erreur : {erreur:.4e}")

    # Sauvegarde des stats globales dans la racine
    fic.create_dataset("L2_error", data=np.array(tab_err))
    fic.create_dataset("CPU_Time", data=np.array(tab_cpu))

print("\nFichier resultats.hdf5 est prêt")

#  Partie Récupération et VTK
def write_vtk(name, N, U_data):
    # On rajoute les bords à 0 pour Paraview
    M = N + 2
    h = 1.0 / (N + 1)
    U_full = np.zeros((M, M))
    U_full[1:-1, 1:-1] = U_data.reshape((N, N))
    
    with open(name, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Resultats Poisson\nASCII\nDATASET STRUCTURED_GRID\n')
        f.write(f'DIMENSIONS {M} {M} 1\n')
        f.write(f'POINTS {M*M} float\n')
        for j in range(M):
            for i in range(M):
                f.write(f"{i*h} {j*h} 0.0\n")
        f.write(f'\nPOINT_DATA {M*M}\nSCALARS solution float 1\nLOOKUP_TABLE default\n')
        for val in U_full.ravel():
            f.write(f"{val}\n")

# Exemple de lecture du HDF5 pour générer un VTK
with h5py.File("../results/resultats.hdf5", "r") as f_in:
    # On récupère les données du groupe N=10
    g10 = f_in["Sol_N10"] 
    n_val = g10.attrs['N']
    u_val = g10["Solution"][()]
    
    write_vtk("../results/visualisation_N10.vtk", n_val, u_val)
    print("Fichier visualisation_N10.vtk prêt pour Paraview.")
