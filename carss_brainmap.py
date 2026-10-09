"""
================================================================
 CARSS  --  Run + Colored Brain Map
================================================================
 Runs CARSS on the two-source data and draws the COLORED cortex:
   - faces colored by the CARSS localizer map mu(v)
   - green dots at found sources
   - red triangles at ground truth
 This is the figure that shows CARSS's result on the brain.
================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.io import loadmat

def load_any(path, keys):
    """Load keys from a .mat file, trying scipy (old format) first,
    then mat73 (v7.3). Handles both without the caller knowing which."""
    out = {}
    try:
        from scipy.io import loadmat as _sl
        m = _sl(path)
        for k in keys:
            out[k] = m[k]
    except (NotImplementedError, ValueError):
        import mat73
        m = mat73.loadmat(path)
        for k in keys:
            out[k] = m[k]
    return out

# ---------- paths ----------
LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_two.mat'
MESHFILE  = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\LFD_Fxd_Vertices.mat'

# ---------- leadfield (v7.3) ----------
K = np.asarray(load_any(LEADFIELD, ['K'])['K'], dtype=float)
n_ch = K.shape[0]; n_dip = K.shape[1]//3
print(f"Leadfield {n_ch} ch, {n_dip} dipoles")

# ---------- data (v7.3) ----------
Dd = load_any(DATA, ['EEG_clean','GT_A_dip','GT_A_pos_mm','GT_B_dip','GT_B_pos_mm','GT_Location'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)
A_dip = int(np.array(Dd['GT_A_dip']).squeeze())-1
B_dip = int(np.array(Dd['GT_B_dip']).squeeze())-1
Location = np.asarray(Dd['GT_Location'], dtype=float)   # 3 x n_dip

# ---------- mesh (OLD format, scipy, squeeze (3,N,1)->(3,N)) ----------
md = loadmat(MESHFILE)
TRI = np.squeeze(np.asarray(md['currytri'])).astype(int)
LOC = np.squeeze(np.asarray(md['curryloc'])).astype(float)
tri = TRI.T - 1
verts = LOC.T
n_faces = tri.shape[0]
print(f"Mesh: {verts.shape[0]} vertices, {n_faces} triangles")

# NOTE: n_faces (29666) should equal n_dip (29666) so per-face = per-dipole
assert n_faces == n_dip, f"faces {n_faces} != dipoles {n_dip}"

# ---------- orthonormal dipole blocks ----------
print("Orthonormalizing dipole blocks...")
Gv_orth = []
for v in range(n_dip):
    q, r = np.linalg.qr(K[:, v*3:v*3+3])
    rank = max(int(np.sum(np.abs(np.diag(r)) > 1e-10)), 1)
    Gv_orth.append(q[:, :rank])

# ---------- CARSS ----------
def subspace_scan(Us, Q=None):
    mu = np.zeros(n_dip); UsT = Us.T
    for v in range(n_dip):
        Go = Gv_orth[v]
        if Q is not None:
            Go = Q @ Go
            q, r = np.linalg.qr(Go)
            rk = int(np.sum(np.abs(np.diag(r)) > 1e-10))
            if rk == 0: mu[v]=0.0; continue
            Go = q[:, :rk]
        s = np.linalg.svd(UsT @ Go, compute_uv=False)
        mu[v] = s[0] if s.size else 0.0
    return mu

C = Y @ Y.T / Y.shape[1]
ev, evec = np.linalg.eigh(C)
ev = ev[::-1]; evec = evec[:, ::-1]; ev = np.maximum(ev,0)
n_eig = int(np.sum(ev/ev[0] > 0.01))
Us = evec[:, :max(n_eig,1)]

# first-iteration localizer map (for coloring the brain) + recursion for peaks
mu_map = subspace_scan(Us, None)     # the localizer landscape to color
found = []; Bcols = np.zeros((n_ch,0))
for it in range(2):                  # two sources
    Q = None if Bcols.shape[1]==0 else (np.eye(n_ch)-Bcols@np.linalg.pinv(Bcols))
    mu = subspace_scan(Us, Q)
    v = int(np.argmax(mu)); found.append(v)
    Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
print(f"Found dipoles: {found}")
print(f"  found[0] err to A = {np.linalg.norm(Location[:,found[0]]-Location[:,A_dip]):.1f} mm")
print(f"  found[1] err to B = {np.linalg.norm(Location[:,found[1]]-Location[:,B_dip]):.1f} mm")

# ---------- draw colored brain ----------
fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection='3d')
polys = verts[tri]

fv = mu_map.copy()
thr = np.percentile(fv, 98)                 # show top 2% in color
norm = (fv - fv.min())/(fv.max()-fv.min()+1e-30)
colors = plt.get_cmap('hot')(norm)
colors[fv < thr] = [0.85,0.85,0.85,0.25]    # faint grey below threshold
pc = Poly3DCollection(polys, facecolors=colors, edgecolors='none')
ax.add_collection3d(pc)

# markers: found (green) + ground truth (red triangles)
for k, d in enumerate(found):
    p = Location[:, d]
    ax.scatter(*p, s=140, c='lime', edgecolors='k', depthshade=False,
               label='found' if k==0 else None)
for k, d in enumerate([A_dip, B_dip]):
    p = Location[:, d]
    ax.scatter(*p, s=200, marker='^', c='red', edgecolors='k', depthshade=False,
               label='ground truth' if k==0 else None)

ax.set_box_aspect([1,1,1]); ax.set_axis_off()
ax.view_init(elev=0, azim=-90)              # adjust to match MATLAB view
ax.set_title('CARSS localizer map (two sources)')
h,l = ax.get_legend_handles_labels()
seen = dict(zip(l,h)); ax.legend(seen.values(), seen.keys(), loc='upper right')
plt.tight_layout()
plt.savefig('carss_result_brain.png', dpi=130)
plt.show()
print("\nSaved: carss_result_brain.png")
print("If the colored hotspots sit on/near the red ground-truth markers,")
print("the brain map is working. Adjust view_init(elev,azim) to match MATLAB.")