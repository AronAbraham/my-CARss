"""
================================================================
 CARSS  --  Step 3b : 3-ORIENTATION Subspace RAP-MUSIC
================================================================
 Faithful version for the UNCONSTRAINED leadfield (3 orientations
 per dipole). Replaces the scalar approximation.

 KEY DIFFERENCE FROM SCALAR:
   Each location v has a 3-column leadfield block Gv (75 x 3),
   spanning a 3D subspace. The localizer measures how well that
   subspace ALIGNS with the signal subspace, via principal angles
   (an SVD between the orthonormalized Gv and the signal projector).

   mu(v) = largest singular value of  (Us' * orth(Gv))
         = cosine of the smallest principal angle between the
           location subspace and the signal subspace.
   A true source -> subspaces align -> mu ~ 1.

 RAP-MUSIC recursion: after finding a source, project its FULL
 3-column block out of the leadfield, then re-scan.

 Expect: iteration 1 -> source A ~0 mm, iteration 2 -> source B,
 now sharper than the scalar version (which put B at 18 mm).
================================================================
"""

import numpy as np

def load_mat(path, keys):
    out = {}
    try:
        from scipy.io import loadmat
        m = loadmat(path)
        for k in keys: out[k] = m[k]
    except NotImplementedError:
        import mat73
        m = mat73.loadmat(path)
        for k in keys: out[k] = np.array(m[k])
    return out

LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_two.mat'

# ---- FULL leadfield (keep all 3 orientations) ----
K = np.asarray(load_mat(LEADFIELD, ['K'])['K'], dtype=float)   # 75 x 88998
n_ch = K.shape[0]
n_dip = K.shape[1] // 3
print(f"Leadfield: {n_ch} channels, {n_dip} dipoles (3 orientations each)")

# ---- two-source data + ground truth ----
Dd = load_mat(DATA, ['EEG_clean','GT_A_dip','GT_A_pos_mm','GT_B_dip','GT_B_pos_mm','GT_Location'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)
A_dip = int(np.array(Dd['GT_A_dip']).squeeze()) - 1
B_dip = int(np.array(Dd['GT_B_dip']).squeeze()) - 1
Location = np.asarray(Dd['GT_Location'], dtype=float)
A_pos = np.array(Dd['GT_A_pos_mm']).squeeze()
B_pos = np.array(Dd['GT_B_pos_mm']).squeeze()
print(f"GT: A dipole {A_dip} at {A_pos}")
print(f"    B dipole {B_dip} at {B_pos}  (sep {np.linalg.norm(A_pos-B_pos):.1f} mm)")

# ---- pre-orthonormalize every dipole's 3-column block ----
# Gv_orth[v] is a 75x3 (or fewer) orthonormal basis for location v's subspace.
print("Orthonormalizing dipole blocks (one-time)...")
Gv_orth = []
for v in range(n_dip):
    Gv = K[:, v*3:v*3+3]                 # 75 x 3
    # economy QR -> orthonormal columns spanning the same subspace
    q, r = np.linalg.qr(Gv)
    # keep only columns with non-negligible norm (rank of the block)
    rank = np.sum(np.abs(np.diag(r)) > 1e-10)
    Gv_orth.append(q[:, :max(rank,1)])
print("  done")

# ================================================================
# 1. Covariance + signal subspace
# ================================================================
C = Y @ Y.T / Y.shape[1]
evals, evecs = np.linalg.eigh(C)
evals = evals[::-1]; evecs = evecs[:, ::-1]
evals = np.maximum(evals, 0)
enorm = evals/evals[0]
n_src = int(np.sum(enorm > 0.01)); n_src = max(n_src, 2)
print(f"\nEigenvalues (top4): {enorm[:4]}")
print(f"Number of sources: {n_src}")
Us = evecs[:, :n_src]                   # 75 x n_src signal subspace

# ================================================================
# 2. Subspace-MUSIC localizer with 3-orientation blocks
#      mu(v) = max singular value of ( Us' * Gv_orth )
#            = cos(smallest principal angle)
#    With optional projector Q applied to the block first (RAP).
# ================================================================
def subspace_scan(Us, Gv_orth, Q=None):
    n_dip = len(Gv_orth)
    mu = np.zeros(n_dip)
    UsT = Us.T
    for v in range(n_dip):
        Go = Gv_orth[v]
        if Q is not None:
            Go = Q @ Go
            # re-orthonormalize after projection
            q, r = np.linalg.qr(Go)
            rank = np.sum(np.abs(np.diag(r)) > 1e-10)
            if rank == 0:
                mu[v] = 0.0; continue
            Go = q[:, :rank]
        # principal angles: singular values of Us' * Go, in [0,1]
        s = np.linalg.svd(UsT @ Go, compute_uv=False)
        mu[v] = s[0] if s.size else 0.0
    return mu

# ================================================================
# 3. RAP-MUSIC recursion (project out FULL blocks)
# ================================================================
found = []
Bcols = np.zeros((n_ch, 0))    # accumulated found-source columns (full blocks)

for it in range(n_src):
    Q = None if Bcols.shape[1] == 0 else (np.eye(n_ch) - Bcols @ np.linalg.pinv(Bcols))
    mu = subspace_scan(Us, Gv_orth, Q)
    v = int(np.argmax(mu))
    found.append(v)
    # add this source's FULL 3-column block to the projector basis
    Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
    p = Location[:, v]
    print(f"\nIter {it+1}: dipole {v}  mu={mu[v]:.4f}  pos [{p[0]:.1f} {p[1]:.1f} {p[2]:.1f}]")
    print(f"   dist to A: {np.linalg.norm(p-A_pos):.1f} mm | to B: {np.linalg.norm(p-B_pos):.1f} mm")

# ================================================================
# 4. Result
# ================================================================
fp = [Location[:, v] for v in found]
errA = min(np.linalg.norm(p-A_pos) for p in fp)
errB = min(np.linalg.norm(p-B_pos) for p in fp)
print(f"\n=== RESULT (3-orientation) ===")
print(f"Best error to A: {errA:.1f} mm")
print(f"Best error to B: {errB:.1f} mm")
if errA < 10 and errB < 10:
    print("STEP 3b PASS -- both sources found with 3-orientation MUSIC.")
    print("Proceed to Step 4 (CARSS certainty layer).")
else:
    print("STEP 3b CHECK -- still missing a source; paste output to debug.")