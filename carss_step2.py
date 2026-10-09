"""
================================================================
 CARSS  --  Step 2 : Basic MUSIC Localizer
================================================================
 This is the REAL validation gate. Step 1 confirmed the subspace
 split; Step 2 confirms it can LOCALIZE.

 It scans every candidate dipole v and computes the MUSIC localizer:

     mu(v) = || Rs * g(v) ||  /  || g(v) ||           (eq. 33)

 where g(v) is that dipole's (scalar) leadfield column and Rs is the
 projector onto the SIGNAL subspace. Sources that lie in the signal
 subspace give mu(v) close to 1; everything else is smaller. The PEAK
 of mu over all dipoles = the estimated source location.

 WHAT TO CHECK:
   The peak dipole should be 29466 (0-based) -- your known source --
   or a dipole essentially at the same position (error ~0 mm).
   If the peak lands on the true source, MUSIC works and CARSS's core
   is validated. We then proceed to Step 3 (RAP-MUSIC recursion for
   multiple sources).

 Run Step 1 first (or this re-does the covariance from scratch).
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

# ---- paths (same as Step 1) ----
LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_aron.mat'

# ---- load leadfield -> scalar G ----
K = np.asarray(load_mat(LEADFIELD, ['K'])['K'], dtype=float)
n_ch, n_dip = K.shape[0], K.shape[1] // 3
G = np.zeros((n_ch, n_dip))
for v in range(n_dip):
    G[:, v] = np.sqrt(np.sum(K[:, v*3:v*3+3]**2, axis=1))

# ---- load data + ground truth ----
Dd = load_mat(DATA, ['EEG_clean', 'GT_ez_dip', 'GT_Location'])
Y  = np.asarray(Dd['EEG_clean'], dtype=float)
gt_dip = int(np.array(Dd['GT_ez_dip']).squeeze()) - 1        # -> 0-based
Location = np.asarray(Dd['GT_Location'], dtype=float)        # 3 x n_dip
gt_pos = Location[:, gt_dip]

# ================================================================
# 1. Covariance + signal subspace (Step 1, condensed)
# ================================================================
C = Y @ Y.T / Y.shape[1]
evals, evecs = np.linalg.eigh(C)
evals = evals[::-1]; evecs = evecs[:, ::-1]
evals = np.maximum(evals, 0)
n_src = int(np.sum(evals/evals[0] > 0.01))
n_src = max(n_src, 1)
Us = evecs[:, :n_src]                       # signal subspace (75 x n_src)
print(f"Signal subspace dimension (n sources): {n_src}")

# ================================================================
# 2. Signal-subspace projector  Rs = Us Us'
# ================================================================
Rs = Us @ Us.T                              # 75 x 75

# ================================================================
# 3. MUSIC localizer over all dipoles
#      mu(v) = ||Rs g(v)|| / ||g(v)||
#    (this is the fraction of g(v)'s energy that lies in the signal
#     subspace; ~1 at a true source, smaller elsewhere)
# ================================================================
print("Scanning all dipoles with MUSIC localizer...")
mu = np.zeros(n_dip)
for v in range(n_dip):
    g = G[:, v]
    ng = np.linalg.norm(g)
    if ng < 1e-30:
        mu[v] = 0.0
    else:
        mu[v] = np.linalg.norm(Rs @ g) / ng

# ================================================================
# 4. Peak = estimated source; compare to ground truth
# ================================================================
peak_dip = int(np.argmax(mu))
peak_pos = Location[:, peak_dip]
err_mm   = np.linalg.norm(peak_pos - gt_pos)

print("\n=== MUSIC RESULT ===")
print(f"Peak dipole (0-based) : {peak_dip}")
print(f"Peak position (mm)    : [{peak_pos[0]:.1f} {peak_pos[1]:.1f} {peak_pos[2]:.1f}]")
print(f"True dipole (0-based) : {gt_dip}")
print(f"True position (mm)    : [{gt_pos[0]:.1f} {gt_pos[1]:.1f} {gt_pos[2]:.1f}]")
print(f"Localization error    : {err_mm:.2f} mm")
print(f"mu at peak            : {mu[peak_dip]:.4f}  (should be ~1)")
print(f"mu at true source     : {mu[gt_dip]:.4f}")

if err_mm < 5:
    print("\nSTEP 2 PASS -- MUSIC localized the source correctly.")
    print("CARSS core is validated. Proceed to Step 3 (RAP-MUSIC recursion).")
else:
    print("\nSTEP 2 CHECK -- peak is not at the true source.")
    print("Do NOT proceed. Paste this output so we debug before building further.")

# ================================================================
# 5. Diagnostic: histogram of mu, and top-10 dipoles
# ================================================================
order = np.argsort(mu)[::-1]
print("\nTop 10 dipoles by localizer value:")
for i in range(10):
    d = order[i]
    p = Location[:, d]
    e = np.linalg.norm(p - gt_pos)
    print(f"  #{i+1}: dipole {d:6d}  mu={mu[d]:.4f}  err_to_truth={e:5.1f} mm")

plt.figure(figsize=(7,4))
plt.hist(mu, bins=100)
plt.axvline(mu[gt_dip], color='r', ls='--', label='mu at true source')
plt.xlabel('MUSIC localizer value mu(v)'); plt.ylabel('# dipoles')
plt.title('Distribution of localizer values\n(true source should be near the max)')
plt.legend(); plt.tight_layout()
plt.savefig('carss_step2_localizer.png', dpi=120)
print("\nSaved: carss_step2_localizer.png")