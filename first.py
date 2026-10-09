"""
================================================================
 CARSS  --  Step 1 : Covariance + Signal/Noise Subspace Split
================================================================
 This is the FIRST validation gate. It does NOT localize yet.
 It only:
   1. loads the real leadfield K (75 x 88998) and collapses it to a
      scalar leadfield G (75 x 29666)  [Option A - validation version]
   2. loads the clean single-source synthetic EEG
   3. computes the measurement covariance C = Y Y'
   4. eigendecomposes C and shows the eigenvalue drop-off
   5. estimates the number of sources n from where eigenvalues plateau

 WHAT TO CHECK:
   For a single noise-free source, the eigenvalue plot should show
   ONE dominant eigenvalue, then a sharp drop to ~zero. That drop
   location = number of sources. If you see one big value then a
   cliff, Step 1 works and we proceed to Step 2 (the MUSIC localizer).

 Requires: numpy, scipy, matplotlib
 Data file: the .mat with EEG_clean, GT_ez_dip, GT_ez_pos_mm, GT_Location
 Leadfield: UnconstCentLFD.mat  (variable K)
================================================================
"""

import numpy as np
import matplotlib.pyplot as plt

# ---- try scipy first; fall back to mat73 for v7.3 files ----
def load_mat(path, keys):
    """Load specified keys from a .mat file, handling both formats."""
    out = {}
    try:
        from scipy.io import loadmat
        m = loadmat(path)
        for k in keys:
            out[k] = m[k]
        print(f"  loaded {path} with scipy.io")
    except NotImplementedError:
        # v7.3 (HDF5) file -> use mat73 or h5py
        try:
            import mat73
            m = mat73.loadmat(path)
            for k in keys:
                out[k] = np.array(m[k])
            print(f"  loaded {path} with mat73 (v7.3)")
        except ImportError:
            raise RuntimeError(
                "This is a v7.3 .mat file. Install mat73:  pip install mat73"
            )
    return out

# ================================================================
# EDIT THESE PATHS
# ================================================================
LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_aron.mat'

# ================================================================
# 1. LOAD LEADFIELD  ->  scalar G (75 x 29666)
# ================================================================
print("Loading leadfield...")
Kd = load_mat(LEADFIELD, ['K'])
K  = np.asarray(Kd['K'], dtype=float)           # 75 x 88998
n_ch  = K.shape[0]
n_dip = K.shape[1] // 3
print(f"  K shape: {K.shape}  ->  {n_ch} channels, {n_dip} dipoles")

# scalar leadfield: norm across the 3 orientation columns per dipole
G = np.zeros((n_ch, n_dip))
for v in range(n_dip):
    cols = slice(v*3, v*3+3)
    G[:, v] = np.sqrt(np.sum(K[:, cols]**2, axis=1))
print(f"  scalar G shape: {G.shape}")

# ================================================================
# 2. LOAD CLEAN SINGLE-SOURCE EEG + GROUND TRUTH
# ================================================================
print("\nLoading synthetic EEG...")
Dd = load_mat(DATA, ['EEG_clean', 'GT_ez_dip', 'GT_ez_pos_mm', 'GT_Location'])
Y  = np.asarray(Dd['EEG_clean'], dtype=float)   # 75 x 300
gt_dip = int(np.array(Dd['GT_ez_dip']).squeeze()) - 1   # MATLAB 1-based -> Python 0-based
gt_pos = np.array(Dd['GT_ez_pos_mm']).squeeze()
Location = np.asarray(Dd['GT_Location'], dtype=float)   # 3 x 29666

print(f"  Y (EEG) shape : {Y.shape}")
print(f"  ground-truth dipole (0-based): {gt_dip}")
print(f"  ground-truth position (mm)   : {gt_pos}")
print(f"  Location shape: {Location.shape}")

# sanity: channel counts must match
assert Y.shape[0] == n_ch, f"EEG has {Y.shape[0]} channels, leadfield has {n_ch}"
assert Location.shape[1] == n_dip, f"Location has {Location.shape[1]} dipoles, leadfield has {n_dip}"

# ================================================================
# 3. MEASUREMENT COVARIANCE  C = Y Y'
# ================================================================
print("\nComputing covariance...")
C = Y @ Y.T / Y.shape[1]        # 75 x 75
print(f"  C shape: {C.shape}")

# ================================================================
# 4. EIGENDECOMPOSITION  ->  signal / noise subspaces
# ================================================================
# C is symmetric PSD; use eigh (returns ascending eigenvalues)
evals, evecs = np.linalg.eigh(C)
evals = evals[::-1]             # descending
evecs = evecs[:, ::-1]
evals = np.maximum(evals, 0)   # clip tiny negatives from numerical error

# normalize for display
evals_norm = evals / evals[0]

print("\nTop 10 eigenvalues (normalized to largest):")
for i in range(min(10, len(evals_norm))):
    print(f"  eig {i+1:2d}: {evals_norm[i]:.4e}")

# ================================================================
# 5. ESTIMATE NUMBER OF SOURCES  n  (where eigenvalues plateau)
# ================================================================
# Simple heuristic: n = number of eigenvalues above a small fraction
# of the largest (e.g. 1%). For a clean single source, expect n = 1.
thresh = 0.01
n_est = int(np.sum(evals_norm > thresh))
print(f"\nEstimated number of sources (eig > {thresh} of max): {n_est}")
print("  (for a clean SINGLE source, this should be 1)")

# signal subspace = top n_est eigenvectors
Us = evecs[:, :max(n_est, 1)]
print(f"  signal subspace Us shape: {Us.shape}")

# ================================================================
# 6. DIAGNOSTIC PLOT  --  the validation gate
# ================================================================
plt.figure(figsize=(8, 5))
plt.subplot(1, 2, 1)
plt.semilogy(range(1, len(evals_norm)+1), evals_norm, 'o-')
plt.axhline(thresh, color='r', ls='--', label=f'{thresh} threshold')
plt.xlabel('Eigenvalue index'); plt.ylabel('Normalized eigenvalue (log)')
plt.title('Eigenvalue drop-off\n(one source -> one big value then cliff)')
plt.legend(); plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(range(1, min(20, len(evals_norm))+1), evals_norm[:20], 'o-')
plt.xlabel('Eigenvalue index'); plt.ylabel('Normalized eigenvalue (linear)')
plt.title('First 20 eigenvalues (linear scale)')
plt.grid(True)
plt.tight_layout()
plt.savefig('carss_step1_eigenvalues.png', dpi=120)
print("\nSaved eigenvalue plot: carss_step1_eigenvalues.png")
print("\n=== STEP 1 CHECK ===")
print("If the eigenvalue plot shows ONE dominant value then a sharp")
print("drop, Step 1 works. Estimated n should be 1 for this clean")
print("single-source case. If so, we proceed to Step 2 (MUSIC localizer).")