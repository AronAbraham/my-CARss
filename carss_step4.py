"""
================================================================
 CARSS  --  Step 4 : Certainty-Based Stopping (the CARSS layer)
================================================================
 Steps 1-3b built RAP-MUSIC. Step 4 adds what makes it CARSS:
 instead of trusting a fixed source count, it RUNS EXTRA iterations,
 records the peak localizer value mu_i at each, and STOPS when mu_i
 drops sharply and plateaus.

 From the paper (eq. 35 + text):
   "The iteration is continued until a significant drop in the
    mu_i(r_i) is observed. After i > n, a similar 'plateau' is
    observed."

 So: real sources -> high mu (~1). Once the true sources are
 exhausted, the next iteration's mu DROPS (no real source left),
 then plateaus low. CARSS keeps sources found BEFORE the drop.

 On the two-source data:
   iter1 mu ~1.0, iter2 mu ~0.84  (two real sources)
   iter3 mu should DROP  (no third source) -> stop, keep 2.

 This DEMONSTRATES the certainty mechanism -- the novel CARSS piece.
================================================================
"""

import numpy as np
import matplotlib.pyplot as plt

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

# ---- full leadfield, orthonormalized blocks ----
K = np.asarray(load_mat(LEADFIELD, ['K'])['K'], dtype=float)
n_ch = K.shape[0]; n_dip = K.shape[1] // 3
print(f"Leadfield: {n_ch} ch, {n_dip} dipoles")

Dd = load_mat(DATA, ['EEG_clean','GT_A_dip','GT_A_pos_mm','GT_B_dip','GT_B_pos_mm','GT_Location'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)
A_pos = np.array(Dd['GT_A_pos_mm']).squeeze()
B_pos = np.array(Dd['GT_B_pos_mm']).squeeze()
Location = np.asarray(Dd['GT_Location'], dtype=float)

print("Orthonormalizing dipole blocks...")
Gv_orth = []
for v in range(n_dip):
    q, r = np.linalg.qr(K[:, v*3:v*3+3])
    rank = max(int(np.sum(np.abs(np.diag(r)) > 1e-10)), 1)
    Gv_orth.append(q[:, :rank])
print("  done")

# ---- covariance + signal subspace ----
C = Y @ Y.T / Y.shape[1]
evals, evecs = np.linalg.eigh(C)
evals = evals[::-1]; evecs = evecs[:, ::-1]; evals = np.maximum(evals, 0)
enorm = evals/evals[0]
n_eig = int(np.sum(enorm > 0.01))
print(f"\nEigenvalue-based source count: {n_eig}")
Us = evecs[:, :max(n_eig,1)]

def subspace_scan(Us, Gv_orth, Q=None):
    mu = np.zeros(len(Gv_orth)); UsT = Us.T
    for v in range(len(Gv_orth)):
        Go = Gv_orth[v]
        if Q is not None:
            Go = Q @ Go
            q, r = np.linalg.qr(Go)
            rank = int(np.sum(np.abs(np.diag(r)) > 1e-10))
            if rank == 0: mu[v] = 0.0; continue
            Go = q[:, :rank]
        s = np.linalg.svd(UsT @ Go, compute_uv=False)
        mu[v] = s[0] if s.size else 0.0
    return mu

# ================================================================
# CERTAINTY-DRIVEN RECURSION
#   Run MAX_ITER (more than n_eig), record peak mu each time,
#   detect the drop.
# ================================================================
MAX_ITER = max(n_eig + 3, 5)     # run a few extra to see the drop
print(f"\nRunning {MAX_ITER} iterations to observe the certainty drop...")

peak_mu = []
found = []
Bcols = np.zeros((n_ch, 0))

for it in range(MAX_ITER):
    Q = None if Bcols.shape[1] == 0 else (np.eye(n_ch) - Bcols @ np.linalg.pinv(Bcols))
    mu = subspace_scan(Us, Gv_orth, Q)
    v = int(np.argmax(mu))
    peak_mu.append(mu[v])
    found.append(v)
    Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
    p = Location[:, v]
    dA = np.linalg.norm(p-A_pos); dB = np.linalg.norm(p-B_pos)
    print(f"  iter {it+1}: dipole {v:6d}  mu={mu[v]:.4f}  dA={dA:5.1f} dB={dB:5.1f}")

# ================================================================
# DETECT THE DROP: largest relative fall between consecutive mu
# ================================================================
peak_mu = np.array(peak_mu)
drops = peak_mu[:-1] - peak_mu[1:]          # positive = mu fell
# certainty count = index of the biggest drop (sources kept = up to there)
biggest_drop_at = int(np.argmax(drops))     # 0-based: drop AFTER this iter
n_certain = biggest_drop_at + 1

print(f"\nPeak mu per iteration: {np.round(peak_mu,4)}")
print(f"Drops between iters  : {np.round(drops,4)}")
print(f"Biggest drop after iteration {n_certain}")
print(f"=> CARSS keeps {n_certain} certain source(s)")

# ---- report the kept sources vs ground truth ----
print("\nKept (certain) sources:")
for i in range(n_certain):
    v = found[i]; p = Location[:, v]
    dA = np.linalg.norm(p-A_pos); dB = np.linalg.norm(p-B_pos)
    print(f"  source {i+1}: dipole {v}  dist A={dA:.1f} mm  B={dB:.1f} mm")

# ---- verdict ----
if n_certain == 2:
    print("\nSTEP 4 PASS -- certainty stopping selected 2 sources (correct).")
    print("The mu drop after iteration 2 is the CARSS certainty signal.")
else:
    print(f"\nSTEP 4 CHECK -- certainty selected {n_certain}, expected 2.")
    print("Look at the mu sequence; the drop may be gentle. Paste output.")

# ================================================================
# PLOT: the certainty curve (this is a figure for your report)
# ================================================================
plt.figure(figsize=(7,4))
plt.plot(range(1, len(peak_mu)+1), peak_mu, 'o-', linewidth=2)
plt.axvline(n_certain+0.5, color='r', ls='--', label=f'stop after {n_certain}')
plt.xlabel('Iteration'); plt.ylabel('Peak localizer value  mu_i')
plt.title('CARSS certainty curve\n(high for real sources, drops+plateaus after)')
plt.legend(); plt.grid(True); plt.tight_layout()
plt.savefig('carss_step4_certainty.png', dpi=120)
print("\nSaved: carss_step4_certainty.png")