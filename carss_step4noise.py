"""
================================================================
 CARSS  --  Step 4 (NOISE), corrected field access
================================================================
 EEG_bySNR loads as a dict with keys snr00/snr03/snr05/snr10.
 Access them directly. Runs the certainty recursion per SNR.
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
        for k in keys: out[k] = m[k]     # keep native (dict for structs)
    return out

LEADFIELD = r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\UnconstCentLFD.mat'
DATA      = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_two.mat'

# ---- leadfield + orthonormal blocks ----
K = np.asarray(load_mat(LEADFIELD, ['K'])['K'], dtype=float)
n_ch = K.shape[0]; n_dip = K.shape[1] // 3
print(f"Leadfield: {n_ch} ch, {n_dip} dipoles")
print("Orthonormalizing dipole blocks (once)...")
Gv_orth = []
for v in range(n_dip):
    q, r = np.linalg.qr(K[:, v*3:v*3+3])
    rank = max(int(np.sum(np.abs(np.diag(r)) > 1e-10)), 1)
    Gv_orth.append(q[:, :rank])
print("  done")

# ---- data ----
Dd = load_mat(DATA, ['EEG_bySNR','GT_A_pos_mm','GT_B_pos_mm','GT_Location'])
A_pos = np.array(Dd['GT_A_pos_mm']).squeeze()
B_pos = np.array(Dd['GT_B_pos_mm']).squeeze()
Location = np.asarray(Dd['GT_Location'], dtype=float)
bysnr = Dd['EEG_bySNR']            # dict: snr00, snr03, snr05, snr10
print("EEG_bySNR keys:", list(bysnr.keys()))

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

def run_carss(Y, max_iter=5):
    Y = np.asarray(Y, dtype=float)
    C = Y @ Y.T / Y.shape[1]
    evals, evecs = np.linalg.eigh(C)
    evals = evals[::-1]; evecs = evecs[:, ::-1]; evals = np.maximum(evals, 0)
    enorm = evals/evals[0]
    n_eig = int(np.sum(enorm > 0.01))
    Us = evecs[:, :max(n_eig,1)]
    peak_mu=[]; found=[]; Bcols=np.zeros((n_ch,0))
    for it in range(max_iter):
        Q = None if Bcols.shape[1]==0 else (np.eye(n_ch) - Bcols @ np.linalg.pinv(Bcols))
        mu = subspace_scan(Us, Gv_orth, Q)
        v = int(np.argmax(mu)); peak_mu.append(mu[v]); found.append(v)
        Bcols = np.column_stack([Bcols, K[:, v*3:v*3+3]])
    peak_mu=np.array(peak_mu); drops=peak_mu[:-1]-peak_mu[1:]
    n_certain=int(np.argmax(drops))+1
    return enorm, n_eig, peak_mu, n_certain, found

# ---- run across SNR (order high->low) ----
order = ['snr10','snr05','snr03','snr00']; labels=[10,5,3,0]
plt.figure(figsize=(8,5)); summary=[]
for name,lab in zip(order,labels):
    if name not in bysnr:
        print(f"[skip] {name} not present"); continue
    Y = np.asarray(bysnr[name], dtype=float)
    enorm,n_eig,peak_mu,n_certain,found = run_carss(Y)
    errs=[(np.linalg.norm(Location[:,found[i]]-A_pos),
           np.linalg.norm(Location[:,found[i]]-B_pos)) for i in range(n_certain)]
    bestA=min(e[0] for e in errs); bestB=min(e[1] for e in errs)
    print(f"\nSNR {lab:2d} dB | eig-count={n_eig} | certainty-count={n_certain}")
    print(f"   peak mu: {np.round(peak_mu,4)}")
    print(f"   best err A={bestA:.1f} mm  B={bestB:.1f} mm")
    summary.append((lab,n_eig,n_certain,bestA,bestB))
    plt.plot(range(1,len(peak_mu)+1), peak_mu,'o-',label=f'{lab} dB')

plt.xlabel('Iteration'); plt.ylabel('Peak localizer mu_i')
plt.title('CARSS certainty curve vs noise')
plt.legend(); plt.grid(True); plt.tight_layout()
plt.savefig('carss_step4_noise.png', dpi=120)
print("\nSaved: carss_step4_noise.png")
print("\n=== SUMMARY ===")
print(" SNR | eig | certain | errA | errB")
for lab,ne,nc,ea,eb in summary:
    print(f" {lab:3d} |  {ne}  |    {nc}    | {ea:4.1f} | {eb:4.1f}")