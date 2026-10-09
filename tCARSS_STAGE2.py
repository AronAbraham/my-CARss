"""
================================================================
 CARSS  --  Stage II : Graph-Regularized Inverse (reduced space)
================================================================
 Completes the two-stage CARSS:
   Stage I  (done): certainty recursion -> the CERTAIN source dipoles
   Stage II (this): estimate source AMPLITUDES on a reduced space
                    (certain dipoles + neighborhoods), with a graph
                    Laplacian smoothness term.

 Reduced-space graph-regularized least squares (closed form):
     j_red = (Kr^T Kr + mu * Lr)^{-1} Kr^T phi
 where Kr, Lr are the leadfield and Laplacian restricted to the
 reduced dipole set, phi is the measurement at the event peak.

 This is faithful to CARSS ("Certainty-based REDUCED Sparse Solution"):
 reduce the space (Stage I), then solve on it (Stage II). The graph
 term (mu*Lr) enforces spatial smoothness within the source patches.

 VALIDATION: estimated amplitudes should peak at the true sources,
 and the map should be spatially smooth (graph term working).
================================================================
"""
import numpy as np
import scipy.sparse as sp
from scipy.io import loadmat

def load_any(path, keys):
    out = {}
    try:
        m = loadmat(path)
        for k in keys: out[k] = m[k]
    except (NotImplementedError, ValueError):
        import mat73
        m = mat73.loadmat(path)
        for k in keys: out[k] = m[k]
    return out

LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_two.mat'
LAPL      = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\Laplacian.mat'
NEIGH     = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\Neighbor.mat'

# ---------- load ----------
K = np.asarray(load_any(LEADFIELD, ['K'])['K'], dtype=float)     # 75 x 88998
n_ch = K.shape[0]; n_dip = K.shape[1]//3
Dd = load_any(DATA, ['EEG_clean','GT_A_dip','GT_A_pos_mm','GT_B_dip','GT_B_pos_mm','GT_Location'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)
A_dip = int(np.array(Dd['GT_A_dip']).squeeze())-1
B_dip = int(np.array(Dd['GT_B_dip']).squeeze())-1
Location = np.asarray(Dd['GT_Location'], dtype=float)
L_full = loadmat(LAPL)['L']
if not sp.issparse(L_full): L_full = sp.csr_matrix(L_full)
L_full = L_full.tocsr().astype(float)          # 29666 x 29666 (per-dipole)
Neigh = loadmat(NEIGH)['Neighbor_inclusive'].astype(int) - 1   # 29666 x 9, 0-based
print(f"K {K.shape}, L {L_full.shape} nnz={L_full.nnz}, Neigh {Neigh.shape}")

# scalar leadfield (per dipole) for the amplitude inverse
G = np.zeros((n_ch, n_dip))
for v in range(n_dip):
    G[:, v] = np.sqrt(np.sum(K[:, v*3:v*3+3]**2, axis=1))

# ---------- Stage I : CARSS certainty recursion (condensed) ----------
# (reuse the validated subspace RAP-MUSIC; find certain dipoles)
print("Orthonormalizing dipole blocks...")
Gv_orth = []
for v in range(n_dip):
    q, r = np.linalg.qr(K[:, v*3:v*3+3])
    rk = max(int(np.sum(np.abs(np.diag(r)) > 1e-10)), 1)
    Gv_orth.append(q[:, :rk])

def scan(Us, Q=None):
    mu = np.zeros(n_dip); UsT = Us.T
    for v in range(n_dip):
        Go = Gv_orth[v]
        if Q is not None:
            Go = Q @ Go
            q, r = np.linalg.qr(Go); rk = int(np.sum(np.abs(np.diag(r))>1e-10))
            if rk == 0: mu[v]=0.0; continue
            Go = q[:, :rk]
        s = np.linalg.svd(UsT@Go, compute_uv=False); mu[v] = s[0] if s.size else 0.0
    return mu

C = Y @ Y.T / Y.shape[1]
ev, evec = np.linalg.eigh(C); ev=ev[::-1]; evec=evec[:, ::-1]; ev=np.maximum(ev,0)
n_eig = int(np.sum(ev/ev[0] > 0.01)); Us = evec[:, :max(n_eig,1)]
found=[]; peak=[]; Bcols=np.zeros((n_ch,0))
for it in range(n_eig+2):
    Q = None if Bcols.shape[1]==0 else (np.eye(n_ch)-Bcols@np.linalg.pinv(Bcols))
    mu = scan(Us, Q); v=int(np.argmax(mu)); peak.append(mu[v]); found.append(v)
    Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
peak=np.array(peak); drops=peak[:-1]-peak[1:]; n_certain=int(np.argmax(drops))+1
certain = found[:n_certain]
print(f"\nStage I: certain dipoles = {certain} (n_certain={n_certain})")

# ---------- build REDUCED space: certain dipoles + neighborhoods ----------
def expand(dips, hops=1):
    s = set(dips)
    for _ in range(hops):
        add = set()
        for d in list(s):
            for nb in Neigh[d]:
                if nb >= 0: add.add(int(nb))
        s |= add
    return sorted(s)

red = expand(certain, hops=2)          # certain sources + 2-hop neighborhood
red = np.array(red)
print(f"Reduced space: {len(red)} dipoles (from {n_certain} certain + neighbors)")

Kr = G[:, red]                          # 75 x |red|  (scalar leadfield reduced)
Lr = L_full[np.ix_(red, red)].toarray() # |red| x |red|  reduced Laplacian (small now)

# ---------- Stage II : graph-regularized inverse ----------
# use the measurement at the peak-power time sample
pk_t = int(np.argmax(np.sum(Y**2, axis=0)))
phi = Y[:, pk_t]                         # 75

mu_graph = 0.1                           # graph weight (tune)
#mu_graph = 0.0                           # graph weight (tune)
#mu_graph = 1.0                           # graph weight (tune)
A = Kr.T @ Kr + mu_graph * Lr
j_red = np.linalg.solve(A, Kr.T @ phi)  # |red| amplitudes

# map back to full brain
j_full = np.zeros(n_dip)
j_full[red] = np.abs(j_red)
j_full /= (j_full.max() + 1e-30)

# ---------- validate ----------
peak_dip = red[int(np.argmax(np.abs(j_red)))]
print(f"\nStage II peak amplitude at dipole {peak_dip}")
print(f"  dist to A: {np.linalg.norm(Location[:,peak_dip]-Location[:,A_dip]):.1f} mm")
print(f"  dist to B: {np.linalg.norm(Location[:,peak_dip]-Location[:,B_dip]):.1f} mm")

# how much amplitude sits near each true source?
for name, td in [('A', A_dip), ('B', B_dip)]:
    d = np.linalg.norm(Location - Location[:, td:td+1], axis=0)
    near = np.where(d < 15)[0]
    amp_near = j_full[near].max() if near.size else 0
    print(f"  max amplitude within 15mm of {name}: {amp_near:.3f}")

# graph smoothness of the estimate (lower = smoother)
smooth = j_red @ (Lr @ j_red)
print(f"\nGraph smoothness j^T L j = {smooth:.4e} (mu_graph={mu_graph})")
print("Compare with mu_graph=0 (no smoothing) to see the graph term's effect.")

np.save('carss_stage2_jfull.npy', j_full)
np.save('carss_stage2_reduced.npy', red)
print("\nSaved source amplitude map -> carss_stage2_jfull.npy")
print("(feed j_full to the brain plotter as face_values to visualize)")

print("\n=== STAGE II CHECK ===")
print("If peak amplitude is near a true source and both A,B have")
print("amplitude nearby, the graph-regularized reduced inverse works.")
print("Try mu_graph = 0, 0.1, 1.0 to see smoothing increase.")