T1 is validated — all three source signatures are visible, and I can point to each one:

- **2-6 Hz peak (the maximum):** source A's **spike** (2 Hz, energy smearing up to 6 Hz) ✓
- **~38-40 Hz peak (the clear second hump):** source B's **40 Hz burst** ✓
- **~110-120 Hz gentle rise:** source A's **120 Hz HFO** — it's there (the little bump around 110-120), just weak ✓

[Certain] All three frequency components your MATLAB generator put in are present in the tensor. T1 correctly transformed the EEG into a space-time-frequency tensor that captures the real spectral content. The wavelet front-end works.

One honest observation about the 120 Hz being weak: [Likely] the HFO shows up small (~0.04 normalized) because in your generator it's scaled at 0.6 and the spike dominates the energy. That's *physically correct* — HFOs are genuinely low-amplitude compared to spikes, which is exactly why they're hard to localize (the whole point of the STSI/HFO literature). So the weak 120 Hz bump isn't a bug; it's the real, expected difficulty of HFOs. Worth noting for your writeup.

[Certain] The `Top-3` printout still says "match ~10-18 Hz?" — ignore that line, it's my stale validation text from the wrong-data assumption. The real validation is what we just confirmed by eye: 2-6 Hz (spike), 40 Hz (source B), 110-120 Hz (HFO) all present. T1 passes.

So the tensor is built and saved (`tcarss_tensor.npy`). **Step T1 done.**

**Next: Step T2 — CPD (Canonical Polyadic Decomposition).** This decomposes your 75×300×70 tensor into R rank-one components, each with a spatial factor (75, over channels), temporal factor (300, over time), and spectral factor (70, over frequency). The validation: the spectral factors should peak at your source frequencies — one component should show a ~40 Hz spectral profile (source B), another the spike's low-freq profile, another the HFO. That's how CPD *separates the sources by their spectral signature* — the core of tCARSS.

Install tensorly first:
```
pip install tensorly
```

Then I'll write Step T2. Before I do, one design question that determines the validation:

[Certain] **How many components (rank R)?** For now, hardcode a reasonable R (say R=4 or 5) to validate that CPD *runs and separates frequencies*. Dynamic rank estimation (auto-choosing R) is Step T4 — we don't need it yet. Starting with fixed R=4 lets us validate the decomposition works before adding the auto-rank complexity. Same discipline: validate the core before the fancy parts.

[Guessing] With R=4 on your data, I'd expect components roughly like: one spike-dominated (low freq), one 40 Hz (source B), one ~120 Hz (HFO), and one noise/residual. If CPD cleanly produces frequency-distinct components, that validates it separates the sources spectrally — which is tCARSS's whole mechanism.

Confirm `pip install tensorly` works, and I'll write Step T2 with fixed R=4, validating that the CPD components separate by frequency. Then T3 (graph regularization), T4 (dynamic rank), and finally wiring the spatial components into your validated CARSS. You're three pieces from complete tCARSS, with the hardest part (CARSS itself) already done.