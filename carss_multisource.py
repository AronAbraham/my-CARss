"""
================================================================
 CARSS  --  Multi-Source Test (synthetic, all Python)
================================================================
 Generates N sources at CONTROLLED minimum separation, all in the
 SAME frequency band but with different time courses (so they are
 distinct, spatially-separated sources -- the paper's test-case
 design). Then runs CARSS and checks it finds all N.

 This is the QUANTITATIVE validation: known ground truth -> real
 localization-error numbers, source-count sweep, noise robustness.

 Design (matches Mannepalli test cases):
   - sources same frequency band, separated only in SPACE
   - controlled minimum inter-source distance
   - the HARD test: method must resolve by spatial subspace alone
     (band-pass could NOT separate these -- only CARSS/MUSIC)

 Uses the real leadfield (UnconstCentLFD.mat). Fully self-contained.
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
        for k in keys: out[k] = m[k]
    return out

LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
MESH      = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\NewMesh.mat'

rng = np.random.default_rng(42)

# ================================================================
# 1. Load leadfield + dipole positions
# ================================================================
K = np.asarray(load_mat(LEADFIELD, ['K'])['K'], dtype=float)   # 75 x 88998
n_ch = K.shape[0]; n_dip = K.shape[1]//3
# dipole positions: New_Mesh(1:3,:) -> 3 x n_dip  (in mm)
mesh = load_mat(MESH, ['New_Mesh'])['New_Mesh']
Location = np.asarray(mesh, dtype=float)[:3, :]    # 3 x n_dip
print(f"Leadfield {n_ch}x{K.shape[1]}, {n_dip} dipoles; positions {Location.shape}")

# pre-orthonormalize dipole blocks (once) -- reused by CARSS
print("Orthonormalizing dipole blocks (once)...")
Gv_orth = []
for v in range(n_dip):
    q, r = np.linalg.qr(K[:, v*3:v*3+3])
    rank = max(int(np.sum(np.abs(np.diag(r)) > 1e-10)), 1)
    Gv_orth.append(q[:, :rank])
print("  done")

# leadfield-strength per dipole (to prefer visible/shallow sources)
Kpow = np.array([np.linalg.norm(K[:, v*3:v*3+3]) for v in range(n_dip)])

# ================================================================
# 2. Choose N well-separated source dipoles
# ================================================================
def choose_sources(N, min_sep_mm, prefer_strong=True):
    order = np.argsort(Kpow)[::-1] if prefer_strong else rng.permutation(n_dip)
    chosen = []
    for v in order:
        p = Location[:, v]
        if all(np.linalg.norm(p - Location[:, c]) >= min_sep_mm for c in chosen):
            chosen.append(int(v))
            if len(chosen) == N:
                break
    if len(chosen) < N:
        raise RuntimeError(f"could only place {len(chosen)} of {N} sources at {min_sep_mm}mm")
    return chosen

# ================================================================
# 3. Build synthetic EEG: N sources, same band, distinct time courses
# ================================================================
def make_eeg(source_dips, sRate=500, n_times=300, patch_mm=8, snr_db=None):
    t = np.arange(n_times)/sRate
    EEG = np.zeros((n_ch, n_times))
    for i, sd in enumerate(source_dips):
        # extended patch around the source dipole
        d = np.linalg.norm(Location - Location[:, sd:sd+1], axis=0)
        patch = np.where(d <= patch_mm)[0]
        taper = np.exp(-0.5*(d[patch]/(patch_mm/2))**2)
        # distinct time course, SAME band (~10-20 Hz), staggered timing
        f0 = 10 + 4*i                 # 10, 14, 18, ... Hz  (all in one band)
        ct = 0.30 + 0.03*i            # slightly staggered
        s = np.exp(-0.5*((t-ct)/0.06)**2) * np.sin(2*np.pi*f0*t)
        for j, v in enumerate(patch):
            EEG += np.outer(K[:, v*3+0], taper[j]*s)   # orientation-1
    if snr_db is not None:
        sig = np.mean(np.var(EEG, axis=1))
        npow = sig/(10**(snr_db/10))
        EEG = EEG + np.sqrt(npow)*rng.standard_normal(EEG.shape)
    return EEG

# ================================================================
# 4. CARSS (subspace RAP-MUSIC + certainty stopping)
# ================================================================
def subspace_scan(Us, Q=None):
    mu = np.zeros(n_dip); UsT = Us.T
    for v in range(n_dip):
        Go = Gv_orth[v]
        if Q is not None:
            Go = Q @ Go
            q, r = np.linalg.qr(Go)
            rank = int(np.sum(np.abs(np.diag(r)) > 1e-10))
            if rank == 0: mu[v]=0.0; continue
            Go = q[:, :rank]
        s = np.linalg.svd(UsT @ Go, compute_uv=False)
        mu[v] = s[0] if s.size else 0.0
    return mu

def carss(Y, max_iter):
    C = Y @ Y.T / Y.shape[1]
    evals, evecs = np.linalg.eigh(C)
    evals = evals[::-1]; evecs = evecs[:, ::-1]; evals = np.maximum(evals,0)
    n_eig = int(np.sum(evals/evals[0] > 0.01))
    Us = evecs[:, :max(n_eig,1)]
    peak=[]; found=[]; Bcols=np.zeros((n_ch,0))
    for it in range(max_iter):
        Q = None if Bcols.shape[1]==0 else (np.eye(n_ch)-Bcols@np.linalg.pinv(Bcols))
        mu = subspace_scan(Us, Q)
        v = int(np.argmax(mu)); peak.append(mu[v]); found.append(v)
        Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
    peak=np.array(peak); drops=peak[:-1]-peak[1:]
    n_certain=int(np.argmax(drops))+1
    return n_eig, peak, n_certain, found

# ================================================================
# 5. Evaluate: for each true source, error to nearest found source
# ================================================================
def eval_result(found, n_certain, true_dips):
    fp = [Location[:, found[i]] for i in range(n_certain)]
    errs = []
    for td in true_dips:
        tp = Location[:, td]
        errs.append(min(np.linalg.norm(tp-p) for p in fp))
    return errs

# ================================================================
# RUN: sweep source count (paper uses 1,3,5)
# ================================================================
print("\n================ MULTI-SOURCE SWEEP ================")
for N, sep in [(1, 0), (3, 40), (5, 40)]:
    true_dips = choose_sources(N, sep) if N > 1 else choose_sources(1, 0)
    print(f"\n--- {N} source(s), min sep {sep} mm ---")
    for snr in [None, 10, 5, 3]:          # clean + noisy
        Y = make_eeg(true_dips, snr_db=snr)
        n_eig, peak, n_certain, found = carss(Y, max_iter=N+3)
        errs = eval_result(found, n_certain, true_dips)
        tag = "clean" if snr is None else f"{snr}dB"
        print(f"  {tag:5s}: eig={n_eig:2d} certain={n_certain} "
              f"| found/true={n_certain}/{N} "
              f"| mean_err={np.mean(errs):.1f}mm max_err={np.max(errs):.1f}mm")

print("\nCheck: does certainty-count match true N? Do errors stay small?")
print("Where it breaks (source count, noise) is your reportable finding.")