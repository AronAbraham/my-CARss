"""
================================================================
 tCARSS  --  Step T1 : Morlet Wavelet -> 3D Tensor
================================================================
 First piece of the tensor front-end. Converts EEG (channels x time)
 into a 3D tensor (channels x time x frequency) using a continuous
 complex Morlet wavelet transform, per channel.

 Matches the paper's eq. 3 (continuous complex Morlet). We take the
 REAL part of the coefficients for the tensor, consistent with the
 tensor-CARSS and STSI papers.

 VALIDATION (the gate):
   Feed synthetic EEG whose sources are at KNOWN frequencies. The
   tensor's frequency dimension should show energy concentrated at
   those frequencies. If a 10-18 Hz source lights up the 10-18 Hz
   slices and not, say, 60 Hz, the wavelet tensor is correct.

 Requires: numpy, pywavelets (pywt), matplotlib
 (tensorly NOT needed yet -- that's Step T2)
================================================================
"""
import numpy as np
import pywt
import matplotlib.pyplot as plt
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

DATA = r'D:\mathlab\github_stsi\STSI\Spikes\Spike_Example_Patient1\Example_Patient1_hfo_spike_data_two.mat'

# ---- load clean two-source EEG ----
Dd = load_any(DATA, ['EEG_clean','GT_sRate'])
Y = np.asarray(Dd['EEG_clean'], dtype=float)     # channels x time
sRate = float(np.array(Dd['GT_sRate']).squeeze())
n_ch, n_time = Y.shape
print(f"EEG: {n_ch} channels x {n_time} samples @ {sRate} Hz")

# ================================================================
# 1. Choose the frequencies to analyze
#    (the paper uses a band; here we scan a reasonable EEG range.
#     Your two-source data has sources ~10-18 Hz, so include that.)
# ================================================================
freqs = np.arange(2, 141, 2.0)      # 4..60 Hz in 2 Hz steps -> 29 freqs
print(f"Analyzing {len(freqs)} frequencies: {freqs[0]}..{freqs[-1]} Hz")

# convert desired frequencies to CWT scales for complex Morlet 'cmor'
# pywt: frequency = central_freq / (scale * sampling_period)
wavelet = 'cmor1.5-1.0'            # complex Morlet (bandwidth-центре)
fc = pywt.central_frequency(wavelet)
dt = 1.0 / sRate
scales = fc / (freqs * dt)

# ================================================================
# 2. Build the tensor: channels x time x frequency
# ================================================================
print("Building wavelet tensor...")
n_freq = len(freqs)
Tensor = np.zeros((n_ch, n_time, n_freq), dtype=float)
for ch in range(n_ch):
    coef, _ = pywt.cwt(Y[ch, :], scales, wavelet, sampling_period=dt)
    # coef shape: (n_freq, n_time), complex. Take real part.
    Tensor[ch, :, :] = np.real(coef).T     # -> (time, freq)
print(f"Tensor shape: {Tensor.shape}  (channels x time x freq)")

# ================================================================
# 3. VALIDATION: where is the energy in frequency?
#    Sum |tensor| over channels and time -> energy per frequency.
#    Peaks should sit at the source frequencies.
# ================================================================
energy_per_freq = np.sum(np.abs(Tensor), axis=(0, 1))     # length n_freq
energy_per_freq /= energy_per_freq.max()

print("\nFrequency energy profile (normalized):")
for f, e in zip(freqs, energy_per_freq):
    bar = '#' * int(e*40)
    print(f"  {f:4.0f} Hz | {e:.3f} {bar}")

peak_freqs = freqs[np.argsort(energy_per_freq)[::-1][:3]]
print(f"\nTop-3 energy frequencies: {peak_freqs} Hz")
print("Check: do these match your source frequencies (~10-18 Hz)?")

# ================================================================
# 4. Plot the frequency energy profile + a channel's time-freq map
# ================================================================
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].plot(freqs, energy_per_freq, 'o-')
ax[0].set_xlabel('Frequency (Hz)'); ax[0].set_ylabel('Normalized energy')
ax[0].set_title('Energy vs frequency\n(peaks should be at source freqs)')
ax[0].grid(True)

# time-frequency map of the strongest channel
ch_energy = np.sum(np.abs(Tensor), axis=(1, 2))
strong_ch = int(np.argmax(ch_energy))
im = ax[1].imshow(np.abs(Tensor[strong_ch]).T, aspect='auto', origin='lower',
                  extent=[0, n_time/sRate, freqs[0], freqs[-1]], cmap='hot')
ax[1].set_xlabel('Time (s)'); ax[1].set_ylabel('Frequency (Hz)')
ax[1].set_title(f'Time-frequency map (channel {strong_ch})')
plt.colorbar(im, ax=ax[1])
plt.tight_layout()
plt.savefig('tcarss_T1_tensor.png', dpi=120)
print("\nSaved: tcarss_T1_tensor.png")

# save the tensor for Step T2 (CPD)
np.save('tcarss_tensor.npy', Tensor)
np.save('tcarss_freqs.npy', freqs)
print("Saved tensor -> tcarss_tensor.npy (for Step T2)")
print("\n=== T1 CHECK ===")
print("If the energy peaks at the source frequencies and the time-freq")
print("map shows activity at the right band/time, the wavelet tensor is")
print("correct. Then install tensorly and we do Step T2 (CPD).")