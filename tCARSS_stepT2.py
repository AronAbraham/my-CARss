"""
================================================================
 tCARSS  --  Step T2 : Canonical Polyadic Decomposition (CPD)
================================================================
 Decomposes the T1 wavelet tensor (channels x time x freq) into R
 rank-one components, each = spatial(75) x temporal(300) x spectral(70).

 This is the piece that SEPARATES sources by their spectral profile.
 Each component's spectral factor should peak at a different source
 frequency -> that's how the tensor front-end pulls the sources apart
 before CARSS localizes them.

 Fixed rank R=4 for now (dynamic rank is Step T4). Expected on your
 data (spike ~2-6 Hz, source B 40 Hz, HFO ~120 Hz, + residual):
   components with spectral peaks near those frequencies.

 VALIDATION: print each component's peak spectral frequency. If the
 components separate into distinct frequency bands matching the
 sources, CPD works and we proceed to T3 (graph regularization).

 Requires: numpy, tensorly, matplotlib
================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
import tensorly as tl
from tensorly.decomposition import parafac

# ---- load the tensor + freq axis from T1 ----
Tensor = np.load('tcarss_tensor.npy')       # (75, 300, 70)
freqs  = np.load('tcarss_freqs.npy')        # (70,)
n_ch, n_time, n_freq = Tensor.shape
print(f"Loaded tensor {Tensor.shape}, freqs {freqs[0]}..{freqs[-1]} Hz")

# ---- normalize tensor (helps CPD conditioning) ----
T = tl.tensor(Tensor / (np.linalg.norm(Tensor) + 1e-30))

# ================================================================
# CPD  (fixed rank R = 4 for validation)
# ================================================================
#R = 4
R=8
print(f"\nRunning CPD (rank {R})...")
# parafac returns weights + factors [spatial(75xR), temporal(300xR), spectral(70xR)]
weights, factors = parafac(T, rank=R, n_iter_max=500, tol=1e-8,
                           init='svd', random_state=0)
S, Ttemp, F = factors        # spatial, temporal, spectral factor matrices
print(f"Factor shapes: spatial {S.shape}, temporal {Ttemp.shape}, spectral {F.shape}")

# reconstruction fit
T_rec = tl.cp_to_tensor((weights, factors))
fit = 1 - np.linalg.norm(T - T_rec) / np.linalg.norm(T)
print(f"CPD reconstruction fit: {fit:.4f}  (1 = perfect)")

# ================================================================
# VALIDATION: each component's dominant spectral frequency
# ================================================================
print("\n--- Component spectral profiles ---")
for r in range(R):
    fp = np.abs(F[:, r])
    fp_norm = fp / (fp.max() + 1e-30)
    peak_f = freqs[np.argmax(fp)]
    # list frequencies where this component has strong energy (>30% of its peak)
    strong = freqs[fp_norm > 0.3]
    print(f"Component {r+1}: peak {peak_f:.0f} Hz | strong band "
          f"{strong.min():.0f}-{strong.max():.0f} Hz")

# ================================================================
# PLOT: spectral factor of each component
# ================================================================
fig, ax = plt.subplots(2, 1, figsize=(9, 8))
for r in range(R):
    fp = np.abs(F[:, r]); fp = fp/(fp.max()+1e-30)
    ax[0].plot(freqs, fp, 'o-', label=f'comp {r+1}')
ax[0].set_xlabel('Frequency (Hz)'); ax[0].set_ylabel('Norm. spectral weight')
ax[0].set_title('CPD component spectral profiles\n(should separate by source frequency)')
ax[0].legend(); ax[0].grid(True)

for r in range(R):
    tp = np.abs(Ttemp[:, r]); tp = tp/(tp.max()+1e-30)
    ax[1].plot(np.arange(n_time)/500.0, tp, label=f'comp {r+1}')
ax[1].set_xlabel('Time (s)'); ax[1].set_ylabel('Norm. temporal weight')
ax[1].set_title('CPD component temporal profiles')
ax[1].legend(); ax[1].grid(True)
plt.tight_layout()
plt.savefig('tcarss_T2_components.png', dpi=120)
print("\nSaved: tcarss_T2_components.png")

# ---- save factors for T3/wiring into CARSS ----
np.save('tcarss_spatial.npy', S)      # 75 x R  (these feed CARSS later)
np.save('tcarss_spectral.npy', F)
np.save('tcarss_temporal.npy', Ttemp)
print("Saved factors -> tcarss_spatial/spectral/temporal.npy")

print("\n=== T2 CHECK ===")
print("If components separate into distinct frequency bands (one near")
print("the spike ~2-6 Hz, one near 40 Hz, one near 120 Hz), CPD is")
print("separating the sources spectrally -> proceed to T3.")
print("If all components look the same / noisy, we adjust rank or init.")