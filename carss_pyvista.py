"""
================================================================
 CARSS  --  Interactive Brain Map (PyVista)
================================================================
 Smooth, rotatable 3D cortex with the CARSS localizer painted on it,
 plus markers for found sources and ground truth.

 REQUIRES:  pip install pyvista
 If PyVista opens a black/empty window, try:  pip install pyvistaqt
 On a normal local Windows Python this should open an interactive
 window you rotate with the mouse (left-drag rotate, scroll zoom).

 Colors each TRIANGLE (= dipole) by the localizer value -- same
 per-face mapping as the MATLAB trisurf figures.
================================================================
"""
import numpy as np
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

# ---------- paths ----------
LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_two.mat'
MESHFILE  = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\LFD_Fxd_Vertices.mat'

# ---------- load ----------
K = np.asarray(load_any(LEADFIELD, ['K'])['K'], dtype=float)
n_ch = K.shape[0]; n_dip = K.shape[1]//3

Dd = load_any(DATA, ['EEG_clean','GT_A_dip','GT_A_pos_mm','GT_B_dip','GT_B_pos_mm','GT_Location'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)
A_dip = int(np.array(Dd['GT_A_dip']).squeeze())-1
B_dip = int(np.array(Dd['GT_B_dip']).squeeze())-1
Location = np.asarray(Dd['GT_Location'], dtype=float)

md = loadmat(MESHFILE)
TRI = np.squeeze(np.asarray(md['currytri'])).astype(int)     # 3 x 29666
LOC = np.squeeze(np.asarray(md['curryloc'])).astype(float)   # 3 x 14775
tri = TRI.T - 1        # (29666,3) 0-based
verts = LOC.T          # (14775,3)
n_faces = tri.shape[0]
print(f"Mesh: {verts.shape[0]} vertices, {n_faces} triangles; dipoles {n_dip}")

# ---------- CARSS (compute localizer map + found sources) ----------
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
            q, r = np.linalg.qr(Go)
            rk = int(np.sum(np.abs(np.diag(r)) > 1e-10))
            if rk == 0: mu[v]=0.0; continue
            Go = q[:, :rk]
        s = np.linalg.svd(UsT @ Go, compute_uv=False)
        mu[v] = s[0] if s.size else 0.0
    return mu

C = Y @ Y.T / Y.shape[1]
ev, evec = np.linalg.eigh(C); ev = ev[::-1]; evec = evec[:, ::-1]; ev = np.maximum(ev,0)
n_eig = int(np.sum(ev/ev[0] > 0.01)); Us = evec[:, :max(n_eig,1)]

mu_map = scan(Us, None)                      # localizer landscape
found = []; Bcols = np.zeros((n_ch,0))
for it in range(2):
    Q = None if Bcols.shape[1]==0 else (np.eye(n_ch)-Bcols@np.linalg.pinv(Bcols))
    mu = scan(Us, Q); v = int(np.argmax(mu)); found.append(v)
    Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
print(f"Found: {found}  errA={np.linalg.norm(Location[:,found[0]]-Location[:,A_dip]):.1f} "
      f"errB={np.linalg.norm(Location[:,found[1]]-Location[:,B_dip]):.1f} mm")

# ================================================================
# PyVista render
# ================================================================
try:
    import pyvista as pv
except ImportError:
    raise SystemExit("PyVista not installed. Run:  pip install pyvista")

# PolyData faces format: [3, i, j, k, 3, i, j, k, ...]
faces = np.hstack([np.full((n_faces,1), 3, dtype=np.int64), tri]).ravel()
mesh = pv.PolyData(verts, faces)

# per-face scalar = localizer value
mesh.cell_data['localizer'] = mu_map

p = pv.Plotter()
p.add_mesh(mesh, scalars='localizer', cmap='hot', smooth_shading=True,
           show_scalar_bar=True, scalar_bar_args={'title': 'CARSS localizer'})

# found sources (green spheres) and ground truth (red spheres)
for d in found:
    p.add_mesh(pv.Sphere(radius=4, center=Location[:, d]), color='lime')
for d in [A_dip, B_dip]:
    p.add_mesh(pv.Sphere(radius=5, center=Location[:, d]), color='red')

p.add_text('CARSS: green=found, red=ground truth', font_size=10)
p.set_background('white')
print("\nOpening interactive window: left-drag rotate, scroll zoom, q to quit.")
p.show()                      # interactive window

# also save a screenshot for the report
p2 = pv.Plotter(off_screen=True)
p2.add_mesh(mesh, scalars='localizer', cmap='hot', smooth_shading=True)
for d in found:
    p2.add_mesh(pv.Sphere(radius=4, center=Location[:, d]), color='lime')
for d in [A_dip, B_dip]:
    p2.add_mesh(pv.Sphere(radius=5, center=Location[:, d]), color='red')
p2.set_background('white')
p2.screenshot('carss_result_pyvista.png')
print("Saved screenshot: carss_result_pyvista.png")