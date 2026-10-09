"""
================================================================
 CARSS  --  Step 3 : RAP-MUSIC Recursion (multiple sources)
================================================================
 Steps 1-2 validated single-source localization. Step 3 adds the
 RECURSION that finds multiple sources one at a time -- the whole
 reason CARSS/RAP-MUSIC exists.

 Mechanism:
   iteration 1: plain MUSIC -> find strongest source, record its
                leadfield topography b1 = g(v1)
   iteration i: project out already-found sources with
                   Q = I - B (B'B)^-1 B'      (B = [b1 ... b_{i-1}])
                and run MUSIC on the projected space:
                   mu_i(v) = ||Rs Q g(v)|| / ||Q g(v)||
                -> find next source
   repeat for n iterations (n = # sources from eigenvalues)

 WHAT TO CHECK (two-source data):
   1. eigenvalues should show n = 2 (two dominant, then drop)
   2. the two recursions should find BOTH sources -- one near A,
      one near B -- each within a few mm.
   If both are found, RAP-MUSIC works and we proceed to Step 4
   (the CARSS certainty layer).
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

# ---- leadfield -> scalar G ----
K = np.asarray(load_mat(LEADFIELD, ['K'])['K'], dtype=float)
n_ch, n_dip = K.shape[0], K.shape[1] // 3
G = np.zeros((n_ch, n_dip))
for v in range(n_dip):
    G[:, v] = np.sqrt(np.sum(K[:, v*3:v*3+3]**2, axis=1))

# ---- two-source data + ground truth ----
Dd = load_mat(DATA, ['EEG_clean','GT_A_dip','GT_A_pos_mm','GT_B_dip','GT_B_pos_mm','GT_Location'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)
A_dip = int(np.array(Dd['GT_A_dip']).squeeze()) - 1     # 0-based
B_dip = int(np.array(Dd['GT_B_dip']).squeeze()) - 1
Location = np.asarray(Dd['GT_Location'], dtype=float)
A_pos = np.array(Dd['GT_A_pos_mm']).squeeze()
B_pos = np.array(Dd['GT_B_pos_mm']).squeeze()

print(f"Ground truth: source A dipole {A_dip} at {A_pos}")
print(f"              source B dipole {B_dip} at {B_pos}")
print(f"              separation {np.linalg.norm(A_pos-B_pos):.1f} mm")

# ================================================================
# 1. Covariance + eigenvalues -> number of sources
# ================================================================
C = Y @ Y.T / Y.shape[1]
evals, evecs = np.linalg.eigh(C)
evals = evals[::-1]; evecs = evecs[:, ::-1]
evals = np.maximum(evals, 0)
enorm = evals / evals[0]

print("\nTop 6 eigenvalues (normalized):")
for i in range(6):
    print(f"  eig {i+1}: {enorm[i]:.4e}")

n_src = int(np.sum(enorm > 0.01))
print(f"\nEstimated number of sources: {n_src}  (expect 2)")
n_src = max(n_src, 2)   # for this test we expect two; guard against under-count
Us = evecs[:, :n_src]
Rs = Us @ Us.T

# ================================================================
# 2. RAP-MUSIC recursion
# ================================================================
def music_scan(Rs, G, Q=None):
    """Return localizer mu(v) over all dipoles, optionally with
    projector Q applied to remove already-found sources."""
    n_dip = G.shape[1]
    mu = np.zeros(n_dip)
    for v in range(n_dip):
        g = G[:, v]
        gg = g if Q is None else Q @ g
        n = np.linalg.norm(gg)
        mu[v] = 0.0 if n < 1e-30 else np.linalg.norm(Rs @ gg) / n
    return mu

found = []           # list of found dipole indices
B = np.zeros((n_ch, 0))   # matrix of found source topographies

for it in range(n_src):
    if B.shape[1] == 0:
        Q = None
    else:
        # projector removing the span of already-found sources
        Q = np.eye(n_ch) - B @ np.linalg.pinv(B)
    mu = music_scan(Rs, G, Q)
    v = int(np.argmax(mu))
    found.append(v)
    B = np.column_stack([B, G[:, v]])      # add this source's topo
    p = Location[:, v]
    dA = np.linalg.norm(p - A_pos)
    dB = np.linalg.norm(p - B_pos)
    print(f"\nIteration {it+1}: found dipole {v}  mu={mu[v]:.4f}")
    print(f"   position [{p[0]:.1f} {p[1]:.1f} {p[2]:.1f}] mm")
    print(f"   distance to A: {dA:.1f} mm | to B: {dB:.1f} mm")

# ================================================================
# 3. Assess: did we find BOTH sources?
# ================================================================
print("\n=== RESULT ===")
found_pos = [Location[:, v] for v in found]
# nearest-truth assignment
errA = min(np.linalg.norm(p - A_pos) for p in found_pos)
errB = min(np.linalg.norm(p - B_pos) for p in found_pos)
print(f"Best error to source A: {errA:.1f} mm")
print(f"Best error to source B: {errB:.1f} mm")

if errA < 10 and errB < 10:
    print("\nSTEP 3 PASS -- RAP-MUSIC found BOTH sources.")
    print("Proceed to Step 4 (CARSS certainty layer).")
else:
    print("\nSTEP 3 CHECK -- one or both sources missed.")
    print("Paste output; we debug the recursion / projector before Step 4.")